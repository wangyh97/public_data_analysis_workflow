/** Fast, bounded streaming Matrix Market reader using Node built-ins only.
 * Full integer/count/order/coordinate validation; no downsampling.
 * All entries contribute to library totals; only predefined genes are retained.
 * Float64 buffers preserve integer counts exactly below Number.MAX_SAFE_INTEGER.
 * Unlike a whole-matrix import, the full 1.31-billion-entry matrix is never stored.
 */
import fs from 'node:fs';import path from 'node:path';import zlib from 'node:zlib';
import {fileURLToPath} from 'node:url';import {once} from 'node:events';
const root=path.resolve(process.argv[2]);
const src=path.join(root,'data','GEO'),out=path.join(root,'derived','panel_v2','GEO_GSE236581');
const genes=JSON.parse(fs.readFileSync(path.join(root,'panel/_shared/extraction_genes.json'),'utf8'));
const lines=name=>zlib.gunzipSync(fs.readFileSync(path.join(src,name))).toString('utf8').trimEnd().split(/\r?\n/);
const features=lines('GSE236581_features.tsv.gz').map(s=>{const a=s.split('\t');return a[1]||a[0];});
const bars=lines('GSE236581_barcodes.tsv.gz').map(s=>s.replace(/^"|"$/g,''));
const present=genes.filter(g=>features.includes(g)),ng=present.length,n=bars.length;
const lookup=new Int32Array(features.length).fill(-1),isMt=new Uint8Array(features.length);
features.forEach((g,i)=>{const j=present.indexOf(g);if(j>=0){if(features.indexOf(g)!==i)throw Error('Duplicate target symbol');lookup[i]=j;}if(g.startsWith('MT-'))isMt[i]=1;});
if(new Set(bars).size!==n)throw Error('Duplicate barcodes');
const E=new Float64Array(n*ng),lib=new Float64Array(n),det=new Uint32Array(n),mt=new Float64Array(n);
let header=true,line='',dims=null,r=0,c=0,value=0,field=0,seen=0,lastC=-1,lastR=-1,nextProgress=100000000;
const started=Date.now();
const input=fs.createReadStream(path.join(src,'GSE236581_counts.mtx.gz')).pipe(zlib.createGunzip({chunkSize:1024*1024}));
for await(const b of input){
 for(let k=0;k<b.length;k++){
  const ch=b[k];
  if(header){
   if(ch===10){if(!line.startsWith('%')){dims=line.trim().split(/\s+/).map(Number);if(dims[0]!==features.length||dims[1]!==n||dims.length!==3)throw Error('Header dimensions mismatch');header=false;console.log('Dimensions',...dims);}line='';}
   else line+=String.fromCharCode(ch);continue;
  }
  if(ch>=48&&ch<=57){value=value*10+ch-48;continue;}
  if(ch===32||ch===9){if(field===0)r=value-1;else if(field===1)c=value-1;else throw Error('Too many fields');field++;value=0;continue;}
  if(ch===13)continue;
  if(ch!==10)throw Error('Non-integer or unexpected character');
  if(field!==2||r<0||r>=features.length||c<0||c>=n||!Number.isSafeInteger(value))throw Error('Invalid coordinate/count');
  if(c<lastC||(c===lastC&&r<=lastR))throw Error('Unsorted or duplicate coordinate');
  lib[c]+=value;if(value>0)det[c]++;if(isMt[r])mt[c]+=value;
  const j=lookup[r];if(j>=0)E[c*ng+j]=value;
  lastC=c;lastR=r;seen++;value=0;field=0;
  if(seen>=nextProgress){console.log('Entries',seen,'/',dims[2]);nextProgress+=100000000;}
 }
}
if(seen!==dims[2]||field!==0)throw Error('Truncated matrix');
const gz=zlib.createGzip({level:1}),dest=fs.createWriteStream(out+'_expression.csv.gz');gz.pipe(dest);
gz.write(['barcode','library','nFeature','percent_mt',...present].join(',')+'\n');
for(let i=0;i<n;i+=1000){
 let rows='';
 for(let j=i;j<Math.min(i+1000,n);j++)rows+=`${bars[j]},${lib[j]},${det[j]},${100*mt[j]/Math.max(lib[j],1)},${E.subarray(j*ng,(j+1)*ng).join(',')}\n`;
 if(!gz.write(rows))await once(gz,'drain');
}
gz.end();await once(dest,'finish');
const result={n_genes:dims[0],n_cells:dims[1],nnz:dims[2],entries_read:seen,gzip_crc_pass:true,integer_nonnegative:true,coordinates_unique_sorted:true,genes_missing:genes.filter(g=>!present.includes(g)),elapsed_seconds:(Date.now()-started)/1000,node_version:process.version};
fs.writeFileSync(path.join(root,'derived','panel_v2','GEO_GSE236581_validation.json'),JSON.stringify(result,null,2));
console.log('DONE',result);
