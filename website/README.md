# HelixForge documentation site

The HelixForge docs, built with [Astro](https://astro.build) +
[Starlight](https://starlight.astro.build). Deploys as a GitHub Pages **project** site at
<https://rcac-bioinformatics.github.io/helixforge/>.

## Develop

```bash
cd website
npm install          # first time only
npm run dev          # http://localhost:4321/helixforge/
```

`npm run dev` and `npm run build` both run `scripts/sync-cli.mjs` first (see below).

## Build & preview

```bash
npm run build        # sync CLI reference + astro build → dist/
npm run preview      # serve the built dist/ locally
```

## Deploy

Pushing to `main` triggers `.github/workflows/deploy-docs.yml`, which builds the site and
publishes `website/dist/` to GitHub Pages. Enable **Settings → Pages → Build and
deployment → Source: GitHub Actions** once, in the repo settings. You can also trigger it
manually from the Actions tab (`workflow_dispatch`).

## Base path (single source of truth)

The deployment base path is a clearly-commented constant at the **top of
[`astro.config.mjs`](./astro.config.mjs)**:

```js
const SITE = 'https://rcac-bioinformatics.github.io';
const BASE = '/helixforge/';
```

The repo has **no `CNAME` / custom domain**, so the site uses the `/helixforge/` project
prefix. If a custom domain is ever added, set `BASE = '/'`, point `SITE` at the domain, and
add the `CNAME` file to `public/`.

### Internal links and the base path

Astro does **not** auto-prepend `base` to absolute links you author. To keep every internal
link working under `/helixforge/`:

- **Markdown/MDX body links** (`[text](/tutorial/star/)`) are rewritten automatically by
  `rehype-base.mjs` (wired into `markdown.rehypePlugins`). Just write them root-absolute.
- **Component props and hero links** (`<LinkCard href>`, `hero.actions[].link`) are rendered
  outside that pipeline, so they must include the base themselves — use
  `` href={`${import.meta.env.BASE_URL}path/`} `` in MDX, or the literal `/helixforge/…` in
  static frontmatter (only the hero in `index.mdx` needs this).

`npm run build` runs `check-links.mjs` as a **postbuild gate**: it fails the build (and the
Actions deploy) if any internal link is missing the base or does not resolve to an emitted
page. Run it alone with `npm run check-links`.

## CLI reference is generated, not hand-written

`scripts/sync-cli.mjs` mirrors the repo's authoritative `../docs/cli/**` tree (produced by
`scripts/gen_cli_reference.py` from the live `click` app) into
`src/content/docs/cli/**`, adding Starlight frontmatter. That directory is **git-ignored** —
never edit it by hand. To update the CLI reference, regenerate `docs/cli/` in the Python
repo and rebuild; the site follows automatically and cannot drift from the code.

## Content layout

```
src/content/docs/
├── index.mdx              Overview (landing / splash)
├── installation.mdx
├── tutorial/              Maize B73 end-to-end walkthrough (8 steps)
├── cli/                   GENERATED from ../docs/cli (git-ignored)
├── benchmarks.mdx
├── limitations.md
└── citing.md
```

Sidebar and site config live in [`astro.config.mjs`](./astro.config.mjs).
