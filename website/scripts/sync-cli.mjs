// ---------------------------------------------------------------------------
// sync-cli.mjs — mirror the authoritative CLI reference into the Starlight site
// ---------------------------------------------------------------------------
// The repo's `docs/cli/**` tree is generated from the live click app by
// `scripts/gen_cli_reference.py` (option tables come straight from each
// command's docstring/epilog). This script COPIES those Markdown files into
// `src/content/docs/cli/**`, adding the Starlight frontmatter Astro needs.
//
// It never edits the source docs, so the site can never drift from the CLI:
// re-run `gen_cli_reference.py`, rebuild the site, and the reference updates.
// Runs automatically before `astro dev` / `astro build` (see package.json).
// ---------------------------------------------------------------------------
import { promises as fs } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const SRC = path.resolve(__dirname, '../../docs/cli'); // repo docs/cli
const DEST = path.resolve(__dirname, '../src/content/docs/cli'); // Starlight

// Human-friendly page titles keyed by the destination slug. Falls back to the
// Markdown H1 (backticks stripped) when a slug isn't listed here.
const TITLES = {
  'index': 'CLI reference',
};

// One-line descriptions surfaced in the sidebar/search where useful.
function frontmatter(title, description) {
  const esc = (s) => String(s).replace(/"/g, '\\"');
  const lines = ['---', `title: "${esc(title)}"`];
  if (description) lines.push(`description: "${esc(description)}"`);
  // These pages are generated — do not offer an "edit this page" link that
  // points at the mirrored copy.
  lines.push('editUrl: false');
  lines.push('---');
  return lines.join('\n') + '\n\n';
}

// Pull the first Markdown H1 and use it as the title; return {title, body}
// with the H1 line removed (Starlight renders the frontmatter title as H1).
function splitTitle(raw) {
  const lines = raw.split('\n');
  let title = null;
  let idx = -1;
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/^#\s+(.*)$/);
    if (m) {
      title = m[1].replace(/`/g, '').trim();
      idx = i;
      break;
    }
  }
  let body = raw;
  if (idx >= 0) {
    body = lines.slice(0, idx).concat(lines.slice(idx + 1)).join('\n');
  }
  return { title: title || 'Reference', body: body.replace(/^\n+/, '') };
}

// The generated docs use a blockquote "Auto-generated …" banner and relative
// links like `reconcile.md` / `utils/align.md`. Rewrite those bare .md links
// to Starlight directory-style routes so internal links resolve on the site.
function rewriteLinks(body) {
  return body.replace(/\]\((?!https?:|\/|#)([^)]+?)\.md(#[^)]*)?\)/g, (full, target, hash) => {
    // target may be like "reconcile", "utils/align", "../reconcile"
    const clean = target.replace(/^\.\//, '');
    return `](${clean}/${hash || ''})`;
  });
}

async function walk(dir) {
  const out = [];
  for (const entry of await fs.readdir(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) out.push(...(await walk(full)));
    else if (entry.name.endsWith('.md')) out.push(full);
  }
  return out;
}

async function main() {
  const files = await walk(SRC);
  let count = 0;
  for (const file of files) {
    const rel = path.relative(SRC, file); // e.g. "utils/align.md" or "README.md"
    // README.md becomes the CLI overview index page.
    const destRel = rel === 'README.md' ? 'index.md' : rel;
    const slug = destRel.replace(/\.md$/, '');
    const raw = await fs.readFile(file, 'utf8');
    const { title: h1Title, body } = splitTitle(raw);
    const title = TITLES[slug] || h1Title;
    const description =
      slug === 'index'
        ? 'Every HelixForge command and subcommand, generated from the live CLI.'
        : undefined;
    const rewritten = rewriteLinks(body);
    const out = frontmatter(title, description) + rewritten.replace(/\n+$/, '') + '\n';
    const destPath = path.join(DEST, destRel);
    await fs.mkdir(path.dirname(destPath), { recursive: true });
    await fs.writeFile(destPath, out, 'utf8');
    count++;
  }
  console.log(`[sync-cli] wrote ${count} CLI reference pages from ${path.relative(process.cwd(), SRC)}`);
}

main().catch((err) => {
  console.error('[sync-cli] failed:', err);
  process.exit(1);
});
