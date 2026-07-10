// ---------------------------------------------------------------------------
// rehype-base — prepend the site base path to root-absolute internal links
// ---------------------------------------------------------------------------
// Astro does NOT rewrite absolute links you author in Markdown/MDX content:
// a link written as `/tutorial/inputs/` ships verbatim, so on a project site
// served under `/helixforge/` it resolves to `/tutorial/inputs/` (404) instead
// of `/helixforge/tutorial/inputs/`. Only Starlight's own chrome (sidebar,
// logo, prev/next) gets the base automatically.
//
// This rehype plugin walks the rendered content HAST and prepends the base to
// every `href`/`src` that is root-absolute (`/…`), not protocol-relative
// (`//…`), and not already base-prefixed. It is idempotent and never touches
// external URLs, anchors, or already-correct links.
//
// NOTE: this runs on Markdown/MDX *content* only. Links passed as component
// props (e.g. <LinkCard href="…">, Starlight hero `actions`) are rendered
// outside this pipeline and must include the base themselves.

const ATTRS = ['href', 'src'];

function normalizeBase(base) {
  // want a leading slash and NO trailing slash, e.g. "/helixforge"
  let b = base || '/';
  if (!b.startsWith('/')) b = '/' + b;
  return b.replace(/\/+$/, '');
}

function walk(node, fn) {
  fn(node);
  if (node.children) for (const child of node.children) walk(child, fn);
}

export default function rehypeBase(options = {}) {
  const base = normalizeBase(options.base);
  const basePrefix = base + '/'; // "/helixforge/"
  return (tree) => {
    walk(tree, (node) => {
      if (node.type !== 'element' || !node.properties) return;
      for (const attr of ATTRS) {
        const value = node.properties[attr];
        if (typeof value !== 'string' || !value) continue;
        if (!value.startsWith('/')) continue; // relative or external
        if (value.startsWith('//')) continue; // protocol-relative
        if (base && (value === base || value.startsWith(basePrefix))) continue; // already based
        node.properties[attr] = base + value;
      }
    });
  };
}
