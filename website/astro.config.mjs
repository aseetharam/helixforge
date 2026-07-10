// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// ---------------------------------------------------------------------------
// DEPLOYMENT BASE PATH
// ---------------------------------------------------------------------------
// HelixForge deploys as a GitHub Pages *project* site for the org repo:
//   https://rcac-bioinformatics.github.io/helixforge/
// so every internal URL must be prefixed with `/helixforge/`. Astro handles
// this automatically once `site` + `base` are set below — use root-relative
// links written WITHOUT the base (Starlight prepends it) or, in Markdown,
// relative links between pages.
//
// No CNAME / custom domain was found in the repo. If one is ever added, set
// BASE = '/' and SITE = 'https://your-domain' and drop the project prefix.
// ---------------------------------------------------------------------------
const SITE = 'https://rcac-bioinformatics.github.io';
const BASE = '/helixforge/';

export default defineConfig({
  site: SITE,
  base: BASE,
  // Disable SmartyPants typographic substitution. On a CLI-heavy docs site it
  // rewrites double-hyphen flags in prose (e.g. `--transl-table`) into an
  // em-dash (`—transl-table`), which is both wrong and reintroduces the very
  // em-dashes we remove elsewhere. Flags must render literally.
  markdown: { smartypants: false },
  // Emit clean directory-style URLs (…/tutorial/reconcile/) that resolve on
  // GitHub Pages without a trailing-slash redirect.
  trailingSlash: 'always',
  integrations: [
    starlight({
      title: 'HelixForge',
      customCss: ['./src/styles/custom.css'],
      description:
        'Evidence-based refinement of deep-learning gene annotations. Helixer is the backbone; RNA-seq and protein evidence refine a model only when they outscore it.',
      tagline: 'Isoform-aware refinement of Helixer annotations via Mikado reconciliation.',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/aseetharam/helixforge',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/aseetharam/helixforge/edit/main/website/',
      },
      lastUpdated: true,
      sidebar: [
        { label: 'Overview', link: '/' },
        { label: 'Installation', link: '/installation/' },
        {
          label: 'Maize tutorial (B73)',
          items: [
            { label: '1. Inputs overview', link: '/tutorial/inputs/' },
            { label: '2. RNA-seq alignment (STAR)', link: '/tutorial/star/' },
            { label: '3. Transcript assembly (StringTie)', link: '/tutorial/stringtie/' },
            { label: '4. Protein evidence (utils)', link: '/tutorial/protein/' },
            { label: '5. Gene prediction (Helixer)', link: '/tutorial/helixer/' },
            { label: '6. Reconcile', link: '/tutorial/reconcile/' },
            { label: '7. Outputs', link: '/tutorial/outputs/' },
            { label: '8. Preflight (doctor)', link: '/tutorial/preflight/' },
          ],
        },
        {
          label: 'CLI reference',
          items: [
            { label: 'Overview', link: '/cli/' },
            { label: 'reconcile', link: '/cli/reconcile/' },
            { label: 'confidence', link: '/cli/confidence/' },
            { label: 'evidence', link: '/cli/evidence/' },
            { label: 'stats', link: '/cli/stats/' },
            { label: 'doctor', link: '/cli/doctor/' },
            { label: 'viz', link: '/cli/viz/' },
            {
              label: 'utils',
              items: [
                { label: 'align', link: '/cli/utils/align/' },
                { label: 'convert', link: '/cli/utils/convert/' },
                { label: 'extract-proteins', link: '/cli/utils/extract-proteins/' },
                { label: 'fetch-db', link: '/cli/utils/fetch-db/' },
                { label: 'filter', link: '/cli/utils/filter/' },
                { label: 'qc', link: '/cli/utils/qc/' },
                { label: 'summarize', link: '/cli/utils/summarize/' },
              ],
            },
            {
              label: 'parallel',
              items: [
                { label: 'plan', link: '/cli/parallel/plan/' },
                { label: 'tasks', link: '/cli/parallel/tasks/' },
                { label: 'aggregate', link: '/cli/parallel/aggregate/' },
                { label: 'suggest', link: '/cli/parallel/suggest/' },
                { label: 'example-sbatch', link: '/cli/parallel/example-sbatch/' },
              ],
            },
          ],
        },
        { label: 'Benchmarks', link: '/benchmarks/' },
        { label: 'Known limitations', link: '/limitations/' },
        { label: 'Citing', link: '/citing/' },
      ],
    }),
  ],
});
