'use client';

import { useCallback, useEffect, useMemo, useState } from 'react';
import katex from 'katex';

type CoefficientRow = {
  index: number;
  line: number;
  sign: '+' | '-' | '0';
  digits: number;
  prefix: string;
  suffix: string;
  sha256: string;
  is_zero: boolean;
};

type CoefficientIndex = {
  schema: string;
  source: string;
  source_bytes: number;
  source_sha256: string;
  coefficient_order: string;
  rows: CoefficientRow[];
};

const DIGITS_PER_BLOCK = 100;
const BLOCKS_PER_PAGE = 40;
const SOURCE_SHA256 =
  '5058792bf79dd594034393954aac53bf801a57f7a36e989c55efefc9c270fd50';

function MathFormula({ latex, display = true }: { latex: string; display?: boolean }) {
  const html = useMemo(
    () =>
      katex.renderToString(latex, {
        displayMode: display,
        throwOnError: false,
        strict: false,
      }),
    [latex, display],
  );

  return (
    <span
      className={display ? 'math-display' : 'math-inline'}
      dangerouslySetInnerHTML={{ __html: html }}
    />
  );
}

function shortHash(value: string) {
  return `${value.slice(0, 12)}…${value.slice(-12)}`;
}

export default function Home() {
  const [indexData, setIndexData] = useState<CoefficientIndex | null>(null);
  const [selectedIndex, setSelectedIndex] = useState(56);
  const [coefficient, setCoefficient] = useState('');
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [page, setPage] = useState(0);
  const [copied, setCopied] = useState('');

  useEffect(() => {
    fetch('/data/f56_Z_index.json')
      .then((response) => {
        if (!response.ok) throw new Error(`Index request failed: ${response.status}`);
        return response.json() as Promise<CoefficientIndex>;
      })
      .then(setIndexData)
      .catch((reason: unknown) =>
        setError(reason instanceof Error ? reason.message : 'Unable to load the index.'),
      );
  }, []);

  useEffect(() => {
    const controller = new AbortController();
    fetch(`/data/coefficients/z${selectedIndex}.txt`, { signal: controller.signal })
      .then((response) => {
        if (!response.ok) throw new Error(`Coefficient request failed: ${response.status}`);
        return response.text();
      })
      .then((value) => {
        setCoefficient(value.trim());
        setLoading(false);
      })
      .catch((reason: unknown) => {
        if (reason instanceof DOMException && reason.name === 'AbortError') return;
        setError(reason instanceof Error ? reason.message : 'Unable to load the coefficient.');
        setLoading(false);
      });
    return () => controller.abort();
  }, [selectedIndex]);

  const row = indexData?.rows[selectedIndex];
  const unsigned = coefficient.startsWith('-') ? coefficient.slice(1) : coefficient;
  const totalBlocks = Math.max(1, Math.ceil(unsigned.length / DIGITS_PER_BLOCK));
  const totalPages = Math.max(1, Math.ceil(totalBlocks / BLOCKS_PER_PAGE));
  const pageStartBlock = page * BLOCKS_PER_PAGE;
  const pageEndBlock = Math.min(totalBlocks, pageStartBlock + BLOCKS_PER_PAGE);
  const visibleBlocks = useMemo(() => {
    const blocks: Array<{ number: number; start: number; end: number; text: string }> = [];
    for (let block = pageStartBlock; block < pageEndBlock; block += 1) {
      const start = block * DIGITS_PER_BLOCK;
      const end = Math.min(unsigned.length, start + DIGITS_PER_BLOCK);
      blocks.push({ number: block + 1, start: start + 1, end, text: unsigned.slice(start, end) });
    }
    return blocks;
  }, [unsigned, pageStartBlock, pageEndBlock]);

  const copy = useCallback(async (text: string, label: string) => {
    await navigator.clipboard.writeText(text);
    setCopied(label);
    window.setTimeout(() => setCopied(''), 1500);
  }, []);

  const choose = (next: number) => {
    const bounded = Math.max(0, Math.min(56, next));
    if (bounded === selectedIndex) return;
    setLoading(true);
    setError('');
    setPage(0);
    setSelectedIndex(bounded);
  };

  return (
    <main className="min-h-screen">
      <header className="site-header">
        <div className="shell header-inner">
          <a className="brand" href="#top" aria-label="Coefficient reader home">
            <span className="brand-mark">SL</span>
            <span>
              <strong>Exact coefficient reader</strong>
              <small>degree 56 · totally real construction</small>
            </span>
          </a>
          <nav aria-label="Source links">
            <a href="https://doi.org/10.5281/zenodo.22317659" target="_blank" rel="noreferrer">
              Zenodo <span className="icon" aria-hidden="true">↗</span>
            </a>
            <a href="https://github.com/jiangshe294-blip/sl2-f13-totally-real" target="_blank" rel="noreferrer">
              <span className="icon" aria-hidden="true">⑂</span> GitHub
            </a>
          </nav>
        </div>
      </header>

      <div className="shell" id="top">
        <section className="intro-grid">
          <div>
            <p className="eyebrow">PUBLIC MATHEMATICAL DATA · JIANG YU</p>
            <h1>The explicit polynomial, made readable.</h1>
            <p className="lede">
              Inspect every exact integer coefficient of the degree-56 polynomial without
              opening a 31.6 MB line of digits. Values are loaded one at a time, divided into
              numbered 100-digit blocks, and accompanied by independent SHA-256 metadata.
            </p>
          </div>
          <div className="formula-card" aria-label="Polynomial definition">
            <MathFormula latex={'Z(X)=\\sum_{i=0}^{56} z_iX^i'} />
            <MathFormula latex={'D=z_{56}>0,\\qquad f(X)=\\frac{Z(X)}{D}\\in\\mathbf{Q}[X]'} />
            <MathFormula latex={'P_{28}(T)=\\sum_{j=0}^{28}\\frac{z_{2j}}{D}T^j,\\qquad f(X)=P_{28}(X^2)'} />
            <details>
              <summary><span className="icon" aria-hidden="true">&lt;/&gt;</span> LaTeX source</summary>
              <pre>{`Z(X)=\\sum_{i=0}^{56} z_iX^i\nD=z_{56}>0,\\qquad f(X)=\\frac{Z(X)}{D}\\in\\mathbf{Q}[X]\nP_{28}(T)=\\sum_{j=0}^{28}\\frac{z_{2j}}{D}T^j,\\qquad f(X)=P_{28}(X^2)`}</pre>
              <a className="tex-download" href="/f56_reader_source.tex" download>
                <span className="icon" aria-hidden="true">↓</span> Download editable .tex source
              </a>
            </details>
          </div>
        </section>

        <section className="reader" aria-labelledby="reader-title">
          <div className="reader-toolbar">
            <div>
              <p className="section-kicker">COEFFICIENT EXPLORER</p>
              <h2 id="reader-title">Choose <MathFormula latex={'z_i'} display={false} /></h2>
            </div>
            <div className="coefficient-nav">
              <button onClick={() => choose(selectedIndex - 1)} disabled={selectedIndex === 0} aria-label="Previous coefficient">
                <span className="nav-glyph" aria-hidden="true">‹</span>
              </button>
              <label>
                <span>Index</span>
                <select value={selectedIndex} onChange={(event) => choose(Number(event.target.value))}>
                  {Array.from({ length: 57 }, (_, index) => (
                    <option key={index} value={index}>z{index}{index % 2 === 1 ? ' = 0' : ''}</option>
                  ))}
                </select>
              </label>
              <button onClick={() => choose(selectedIndex + 1)} disabled={selectedIndex === 56} aria-label="Next coefficient">
                <span className="nav-glyph" aria-hidden="true">›</span>
              </button>
            </div>
          </div>

          {error ? <div className="error-box">{error}</div> : null}

          <div className="metadata-grid">
            <div><span>Coefficient</span><strong>z<sub>{selectedIndex}</sub></strong></div>
            <div><span>Source line</span><strong>{row?.line ?? '—'}</strong></div>
            <div><span>Sign</span><strong>{row?.sign ?? '—'}</strong></div>
            <div><span>Decimal digits</span><strong>{row?.digits.toLocaleString() ?? '—'}</strong></div>
            <div className="hash-cell"><span><span className="icon" aria-hidden="true">#</span> SHA-256</span><code title={row?.sha256}>{row ? shortHash(row.sha256) : '—'}</code></div>
          </div>

          <div className="value-actions">
            <div>
              <p className="value-label">Exact integer value</p>
              <p className="value-status">
                {loading ? 'Loading exact digits…' : row?.is_zero ? 'The coefficient is exactly zero.' : `${totalBlocks.toLocaleString()} numbered blocks`}
              </p>
            </div>
            <div className="action-group">
              <button className="action-button" disabled={loading || !coefficient} onClick={() => copy(coefficient, 'integer')}>
                <span className="icon" aria-hidden="true">{copied === 'integer' ? '✓' : '⧉'}</span>
                {copied === 'integer' ? 'Copied' : 'Copy integer'}
              </button>
              <button className="action-button" disabled={loading || !coefficient} onClick={() => copy(`(${coefficient})/D`, 'rational')}>
                <span className="icon" aria-hidden="true">{copied === 'rational' ? '✓' : '⧉'}</span>
                {copied === 'rational' ? 'Copied' : 'Copy zᵢ/D'}
              </button>
              <a className="action-button" href={`/data/coefficients/z${selectedIndex}.txt`} download>
                <span className="icon" aria-hidden="true">↓</span> Download .txt
              </a>
            </div>
          </div>

          <div className="digit-panel" aria-busy={loading}>
            {loading ? (
              <div className="loading-line"><span />Loading coefficient z{selectedIndex}…</div>
            ) : row?.is_zero ? (
              <div className="zero-value">0</div>
            ) : (
              <>
                {coefficient.startsWith('-') ? <div className="minus-sign">−</div> : null}
                <div className="blocks">
                  {visibleBlocks.map((block) => (
                    <div className="digit-row" key={block.number}>
                      <span className="block-number">{block.number.toLocaleString()}</span>
                      <code>{block.text}</code>
                      <span className="block-offset">digits {block.start.toLocaleString()}–{block.end.toLocaleString()}</span>
                    </div>
                  ))}
                </div>
              </>
            )}
          </div>

          {!loading && !row?.is_zero ? (
            <div className="page-nav">
              <button onClick={() => setPage((value) => Math.max(0, value - 1))} disabled={page === 0}>
                <span className="icon" aria-hidden="true">←</span> Previous blocks
              </button>
              <label>
                Page
                <input
                  type="number"
                  min={1}
                  max={totalPages}
                  value={page + 1}
                  onChange={(event) => {
                    const value = Number(event.target.value);
                    if (Number.isFinite(value)) setPage(Math.max(0, Math.min(totalPages - 1, value - 1)));
                  }}
                />
                of {totalPages.toLocaleString()}
              </label>
              <button onClick={() => setPage((value) => Math.min(totalPages - 1, value + 1))} disabled={page === totalPages - 1}>
                Next blocks <span className="icon" aria-hidden="true">→</span>
              </button>
            </div>
          ) : null}
        </section>

        <section className="integrity" aria-labelledby="integrity-title">
          <div>
            <p className="section-kicker">REPRODUCIBILITY</p>
            <h2 id="integrity-title">A readable view, not a replacement dataset.</h2>
            <p>
              The canonical file remains <code>f56_Z.txt</code>. Line <MathFormula latex={'i+1'} display={false} /> is exactly
              the coefficient <MathFormula latex={'z_i'} display={false} />, in ascending order. This site splits those same
              lines into separate files only to make browser access practical.
            </p>
          </div>
          <dl>
            <div><dt>Canonical size</dt><dd>31,560,750 bytes</dd></div>
            <div><dt>Canonical SHA-256</dt><dd><code>{SOURCE_SHA256}</code></dd></div>
            <div><dt>Order</dt><dd>ascending, z<sub>0</sub> through z<sub>56</sub></dd></div>
            <div><dt>License</dt><dd>CC BY 4.0</dd></div>
          </dl>
        </section>

        <footer>
          <p>Exact coefficient companion to the construction by Jiang Yu.</p>
          <div>
            <a href="https://doi.org/10.5281/zenodo.22317659" target="_blank" rel="noreferrer">DOI 10.5281/zenodo.22317659</a>
            <span>·</span>
            <a href="https://github.com/jiangshe294-blip/sl2-f13-totally-real" target="_blank" rel="noreferrer">Source and verification</a>
          </div>
        </footer>
      </div>
    </main>
  );
}
