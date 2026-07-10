// Verify internal links: every root-absolute internal link must (a) carry the
// base prefix and (b) resolve to an emitted file in dist/.
import { promises as fs } from 'node:fs';
import path from 'node:path';
const DIST = 'dist', BASE = '/helixforge/';
async function walk(d){const o=[];for(const e of await fs.readdir(d,{withFileTypes:true})){const f=path.join(d,e.name);if(e.isDirectory())o.push(...await walk(f));else if(e.name.endsWith('.html'))o.push(f);}return o;}
const exists=async p=>{try{await fs.access(p);return true;}catch{return false;}};
const files=await walk(DIST), missingBase=[], broken=[];
for(const file of files){
  const html=await fs.readFile(file,'utf8');
  for(const m of html.matchAll(/(?:href|src)="([^"]+)"/g)){
    let href=m[1];
    if(!href.startsWith('/')||href.startsWith('//')) continue; // external/relative
    if(!href.startsWith(BASE)){ missingBase.push(`${path.relative(DIST,file)}  ->  ${href}`); continue; }
    href=href.split('#')[0].split('?')[0];
    if(href.slice(BASE.length).startsWith('pagefind')||href.slice(BASE.length).startsWith('_astro')) continue;
    const rel=href.slice(BASE.length);
    const target=href.endsWith('/')?path.join(DIST,rel,'index.html'):path.join(DIST,rel);
    if(!(await exists(target))&&!(await exists(path.join(DIST,rel,'index.html')))) broken.push(`${path.relative(DIST,file)}  ->  ${m[1]}`);
  }
}
const u=a=>[...new Set(a)];
if(missingBase.length){console.log('MISSING BASE PREFIX:\n'+u(missingBase).join('\n'));}
if(broken.length){console.log('BROKEN (no target file):\n'+u(broken).join('\n'));}
if(!missingBase.length&&!broken.length) console.log(`OK: all internal links carry the base and resolve, across ${files.length} pages`);
else process.exit(1);
