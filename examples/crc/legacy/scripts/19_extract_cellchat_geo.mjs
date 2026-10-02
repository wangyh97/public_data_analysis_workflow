/** Full-gene sparse extraction for the frozen CellChat cells; uint32 little-endian
 * triplets (gene index, local cell index, count). No dense full matrix is created.
 * The complete original matrix is streamed and its coordinate order is checked.
 */
import fs from 'node:fs';import path from 'node:path';import zlib from 'node:zlib';import {fileURLToPath} from 'node:url';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..'),out=path.join(root,'data/CellChat');
const csv=fs.readFileSync(path.join(out,'selected_cells.csv'),'utf8').trim().split(/\r?\n/).map(x=>x.split(','));const names=csv.shift();
const rows=csv.map(x=>Object.fromEntries(names.map((n,i)=>[n,x[i]]))).filter(x=>x.dataset==='GEO');
const mapping=new Map(rows.map(x=>[+x.matrix_col,{job:x.job,col:+x.local_col}]));
const jobs=[...new Set(rows.map(x=>x.job))],state={};
const genes=zlib.gunzipSync(fs.readFileSync(path.join(root,'data/GEO/GSE236581_features.tsv.gz'))).toString().trimEnd().split(/\r?\n/).map(x=>x.split('\t')[1]);
for(const j of jobs){state[j]={fd:fs.openSync(path.join(out,j+'_counts.bin'),'w'),buf:Buffer.alloc(12*65536),used:0,nnz:0};fs.writeFileSync(path.join(out,j+'_genes.txt'),genes.join('\n')+'\n');}
function emit(j,r,c,v){const s=state[j];s.buf.writeUInt32LE(r,s.used);s.buf.writeUInt32LE(c,s.used+4);s.buf.writeUInt32LE(v,s.used+8);s.used+=12;s.nnz++;if(s.used===s.buf.length){fs.writeSync(s.fd,s.buf);s.used=0;}}
let header=true,line='',dims,r=0,c=0,v=0,field=0,seen=0,lastC=0,lastR=0;const start=Date.now();
for await(const b of fs.createReadStream(path.join(root,'data/GEO/GSE236581_counts.mtx.gz')).pipe(zlib.createGunzip({chunkSize:1048576}))){
 for(const ch of b){
  if(header){if(ch===10){if(!line.startsWith('%')){dims=line.trim().split(/\s+/).map(Number);header=false;}line='';}else line+=String.fromCharCode(ch);continue;}
  if(ch>=48&&ch<=57){v=v*10+ch-48;continue;}
  if(ch===32||ch===9){if(field===0)r=v;else if(field===1)c=v;else throw Error('Invalid fields');field++;v=0;continue;}
  if(ch===13)continue;
  if(ch!==10||field!==2||r<1||r>dims[0]||c<1||c>dims[1]||v>4294967295||c<lastC||(c===lastC&&r<=lastR))throw Error('Invalid coordinate');
  const selected=mapping.get(c);if(selected)emit(selected.job,r,selected.col,v);
  lastC=c;lastR=r;v=0;field=0;seen++;if(seen%200000000===0)console.log('Read',seen,'/',dims[2]);
 }
}
if(seen!==dims[2]||field)throw Error('Truncated');
for(const [job,s]of Object.entries(state)){if(s.used)fs.writeSync(s.fd,s.buf.subarray(0,s.used));fs.closeSync(s.fd);fs.writeFileSync(path.join(out,job+'_extraction.json'),JSON.stringify({job,genes:genes.length,nnz:s.nnz,cells:rows.filter(x=>x.job===job).length,full_matrix_entries_validated:seen,elapsed_seconds:(Date.now()-start)/1000}));}
console.log('COMPLETE',jobs.length,'full-gene sample matrices');
