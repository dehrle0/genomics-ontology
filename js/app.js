/**
 * Genomic Ontology Explorer — app logic
 * Vanilla JS, no build step, no framework. Reads from window.GENES /
 * window.ONTOLOGIES / window.PRS / window.PGX / window.REPORT / window.JOB_META
 * (see data/mock-data.js). Swap that file for real API calls and this file
 * needs no structural changes — see the seam comment at the top of it.
 *
 * Revision 2 additions: Graph layout mode for the ontology tree, a
 * Studies tab (aggregated evidence, separate from curated Publications),
 * reference-links row on Overview, extra variant score columns (CADD,
 * SpliceAI, AlphaMissense, QUAL, read support), and a generated DNA
 * watermark.
 */

(function () {
  "use strict";

  // -----------------------------------------------------------------
  // Shared state
  // -----------------------------------------------------------------
  const state = {
    tree: "hpo",
    layout: "list",              // default to List per user preference
    scopeFindingsOnly: false,
    treeSearch: "",
    selectedGene: null,
    expandedNodes: new Set(),      // ontology group/term ids currently expanded (tree mode)
    expandedVariantRows: new Set() // "gene:variantId" currently expanded (study detail)
  };

  const byId = (id) => document.getElementById(id);
  const geneBySymbol = (sym) => GENES.find((g) => g.symbol === sym);
  const esc = (s) => (s === null || s === undefined ? "" : String(s));

  // -----------------------------------------------------------------
  // Top-level tab routing
  // -----------------------------------------------------------------
  function initTabs() {
    document.querySelectorAll(".tabs button").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.querySelectorAll(".tabs button").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        document.querySelectorAll(".view").forEach((v) => v.classList.remove("active"));
        byId("view-" + btn.dataset.view).classList.add("active");
      });
    });
    byId("job-meta").textContent =
      JOB_META.sample + " · OpenCRAVAT " + JOB_META.opencravatVersion + " · " +
      JOB_META.uniqueVariants.toLocaleString() + " variants";
  }

  // -----------------------------------------------------------------
  // ONTOLOGY VIEW — left panel controls
  // -----------------------------------------------------------------
  function initOntologySwitch() {
    document.querySelectorAll("#ontology-switch button").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.querySelectorAll("#ontology-switch button").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        state.tree = btn.dataset.tree;
        state.expandedNodes.clear();
        renderLeftPanel();
      });
    });

    document.querySelectorAll('input[name="layout"]').forEach((r) => {
      r.addEventListener("change", (e) => {
        state.layout = e.target.value;
        renderLeftPanel();
      });
    });

    byId("scope-findings-only").addEventListener("change", (e) => {
      state.scopeFindingsOnly = e.target.checked;
      renderLeftPanel();
    });

    byId("tree-search").addEventListener("input", (e) => {
      state.treeSearch = e.target.value.trim().toLowerCase();
      renderLeftPanel();
    });
  }

  function matchesSearch(text) {
    if (!state.treeSearch) return true;
    return (text || "").toLowerCase().includes(state.treeSearch);
  }

  // A group/term "has findings" if any of its genes have a variant tagged
  // concern / protective / uncertain (i.e. something evaluated, not just
  // "research available").
  function geneHasFindings(symbol) {
    const g = geneBySymbol(symbol);
    if (!g) return false;
    return g.variants.some((v) => v.category !== "uncategorized");
  }

  // Returns the filtered ontology structure shared by Tree/List/Graph:
  // [{ group, terms:[{ term, genes:[...] }] }]  — genes/terms already
  // filtered by search + scope so every renderer draws the same set.
  function filteredOntology() {
    const ont = ONTOLOGIES[state.tree];
    const out = [];
    ont.groups.forEach((group) => {
      const allGroupGenes = Array.from(new Set(group.genes));
      const groupHasFindings = allGroupGenes.some(geneHasFindings);
      if (state.scopeFindingsOnly && !groupHasFindings) return;

      const terms = (group.terms || [])
        .map((term) => {
          let genes = term.genes.slice();
          if (state.scopeFindingsOnly) genes = genes.filter(geneHasFindings);
          if (state.treeSearch) genes = genes.filter(matchesSearch);
          const termHit = matchesSearch(term.label) || matchesSearch(term.id) || genes.length > 0 || term.genes.some(matchesSearch);
          if (state.treeSearch && !termHit) return null;
          if (!genes.length && state.treeSearch) {
            // term text matched but no individual gene matched search -> keep original gene set
            genes = term.genes.slice();
          }
          return { term, genes };
        })
        .filter(Boolean);

      const groupHit =
        matchesSearch(group.label) || matchesSearch(group.id) ||
        allGroupGenes.some(matchesSearch) || terms.length > 0;
      if (state.treeSearch && !groupHit) return;
      if (state.scopeFindingsOnly && terms.every((t) => !t.genes.length) && !groupHasFindings) return;

      out.push({ group, terms });
    });
    return out;
  }

  function renderLeftPanel() {
    if (state.layout === "graph") {
      byId("tree-scroll").style.display = "none";
      renderGraph();
    } else {
      let graphHost = byId("graph-host");
      if (graphHost) graphHost.remove();
      byId("tree-scroll").style.display = "";
      renderTreeOrList();
    }
  }

  // -----------------------------------------------------------------
  // Tree / List rendering (shares one filtered structure; List just
  // skips the intermediate term layer, matching what the user liked)
  // -----------------------------------------------------------------
  function renderTreeOrList() {
    const data = filteredOntology();
    const scroll = byId("tree-scroll");
    scroll.className = "tree-scroll" + (state.layout === "list" ? " tree-list-mode" : "");
    scroll.innerHTML = "";

    data.forEach(({ group, terms }) => {
      const totalGenes = Array.from(new Set(group.genes)).length;
      const groupEl = document.createElement("div");
      groupEl.className = "tree-group";

      const isOpen = state.expandedNodes.has(group.id) || !!state.treeSearch;
      const row = document.createElement("div");
      row.className = "tree-row tree-row--organ" + (isOpen ? " expanded" : "");
      row.innerHTML =
        '<span class="caret">\u25B6</span><span class="node-dot"></span>' +
        "<span>" + group.label + "</span>" +
        '<span class="count-chip">' + totalGenes + " gene" + (totalGenes === 1 ? "" : "s") + "</span>";
      row.addEventListener("click", () => {
        if (state.expandedNodes.has(group.id)) state.expandedNodes.delete(group.id);
        else state.expandedNodes.add(group.id);
        renderTreeOrList();
      });
      groupEl.appendChild(row);

      const childWrap = document.createElement("div");
      childWrap.className = "tree-children" + (isOpen ? " open" : "");

      if (state.layout === "tree" && terms.length) {
        terms.forEach(({ term, genes }) => {
          if (!genes.length) return;
          const termOpen = state.expandedNodes.has(term.id) || !!state.treeSearch;
          const termRow = document.createElement("div");
          termRow.className = "tree-row tree-row--term" + (termOpen ? " expanded" : "");
          termRow.innerHTML =
            '<span class="caret">\u25B6</span><span class="node-dot"></span>' +
            "<span>" + term.label + '</span><span class="count-chip mono">' + term.id + "</span>";
          termRow.addEventListener("click", (e) => {
            e.stopPropagation();
            if (state.expandedNodes.has(term.id)) state.expandedNodes.delete(term.id);
            else state.expandedNodes.add(term.id);
            renderTreeOrList();
          });
          childWrap.appendChild(termRow);

          const geneWrap = document.createElement("div");
          geneWrap.className = "tree-children" + (termOpen ? " open" : "");
          genes.forEach((sym) => geneWrap.appendChild(makeGeneRow(sym)));
          childWrap.appendChild(geneWrap);
        });
      } else {
        // List mode: organ -> genes directly (dedup across terms)
        const flatGenes = Array.from(new Set(
          terms.length ? terms.flatMap((t) => t.genes) : group.genes.filter((s) => !state.treeSearch || matchesSearch(s))
        ));
        flatGenes.forEach((sym) => childWrap.appendChild(makeGeneRow(sym)));
      }

      groupEl.appendChild(childWrap);
      scroll.appendChild(groupEl);
    });

    if (!scroll.children.length) {
      scroll.innerHTML = '<div style="padding:20px;color:var(--slate);font-size:12.5px;">No matches under the current filters.</div>';
    }
  }

  function makeGeneRow(symbol) {
    const g = geneBySymbol(symbol);
    const row = document.createElement("div");
    row.className = "tree-row tree-row--gene" + (state.selectedGene === symbol ? " selected" : "");
    const flagged = g && g.variants.some((v) => v.category === "concern");
    row.innerHTML =
      '<span class="node-dot"></span><span>' + symbol + "</span>" +
      (flagged ? '<span class="count-chip" style="color:var(--concern);border-color:var(--concern-bg);">flag</span>' : "");
    row.addEventListener("click", (e) => {
      e.stopPropagation();
      selectGene(symbol);
    });
    return row;
  }

  function selectGene(symbol) {
    state.selectedGene = symbol;
    if (state.layout === "graph") renderGraph(); else renderTreeOrList();
    renderGeneDetail();
  }

  // -----------------------------------------------------------------
  // Graph layout mode — 2-level node/edge diagram (Group -> Term),
  // with gene chips attached at each term row. Keeps gene selection
  // interactive without the node-dedup complexity of a full 3-tier
  // force graph.
  // -----------------------------------------------------------------
  function renderGraph() {
    const panel = document.querySelector(".tree-panel");
    let host = byId("graph-host");
    if (!host) {
      host = document.createElement("div");
      host.id = "graph-host";
      host.className = "graph-wrap";
      panel.appendChild(host);
    }
    host.innerHTML = "";

    const data = filteredOntology();
    const ROW_H = 28, GROUP_X = 26, TERM_X = 160, GENE_X = 300, PAD_TOP = 20, WIDTH = 620;

    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    svg.setAttribute("class", "graph-svg");
    svg.setAttribute("width", WIDTH);

    const geneChips = []; // { y, genes:[] }
    let y = PAD_TOP;
    const edges = [];
    const nodes = [];

    data.forEach(({ group, terms }) => {
      const groupY = y;
      nodes.push({ type: "group", x: GROUP_X, y: groupY, label: group.label });
      y += ROW_H;

      const visibleTerms = terms.filter((t) => t.genes.length);
      visibleTerms.forEach(({ term, genes }) => {
        const termY = y;
        nodes.push({ type: "term", x: TERM_X, y: termY, label: term.label });
        edges.push({ x1: GROUP_X + 4, y1: groupY, x2: TERM_X - 4, y2: termY });
        geneChips.push({ y: termY, genes });
        y += ROW_H;
      });
      y += 6; // gap between groups
    });

    const totalHeight = y + 10;
    svg.setAttribute("height", totalHeight);
    svg.setAttribute("viewBox", "0 0 " + WIDTH + " " + totalHeight);

    edges.forEach((e) => {
      const line = document.createElementNS(svgns, "line");
      line.setAttribute("x1", e.x1); line.setAttribute("y1", e.y1);
      line.setAttribute("x2", e.x2); line.setAttribute("y2", e.y2);
      line.setAttribute("class", "graph-edge");
      svg.appendChild(line);
    });

    nodes.forEach((n) => {
      const circle = document.createElementNS(svgns, "circle");
      circle.setAttribute("cx", n.x); circle.setAttribute("cy", n.y);
      circle.setAttribute("r", n.type === "group" ? 5 : 4);
      circle.setAttribute("fill", n.type === "group" ? "#1c5b5e" : "#2a5c99");
      svg.appendChild(circle);

      const text = document.createElementNS(svgns, "text");
      text.setAttribute("x", n.x + 10);
      text.setAttribute("y", n.y + 4);
      text.setAttribute("class", "graph-node-label " + n.type);
      text.textContent = n.label;
      svg.appendChild(text);
    });

    host.appendChild(svg);

    const geneLayer = document.createElement("div");
    geneLayer.className = "graph-genes-layer";
    geneLayer.style.width = WIDTH + "px";
    geneLayer.style.height = totalHeight + "px";
    geneChips.forEach((row) => {
      row.genes.forEach((sym, i) => {
        const chip = document.createElement("div");
        chip.className = "graph-gene-chip" + (state.selectedGene === sym ? " selected" : "");
        chip.textContent = sym;
        chip.style.left = (GENE_X + i * 66) + "px";
        chip.style.top = (row.y - 9) + "px";
        chip.addEventListener("click", () => selectGene(sym));
        geneLayer.appendChild(chip);
      });
    });
    host.appendChild(geneLayer);

    if (!data.length) {
      host.innerHTML = '<div style="padding:20px;color:var(--slate);font-size:12.5px;">No matches under the current filters.</div>';
    }
  }

  // -----------------------------------------------------------------
  // ONTOLOGY VIEW — right gene detail
  // -----------------------------------------------------------------
  function renderGeneDetail() {
    const wrap = byId("gene-detail");
    const g = geneBySymbol(state.selectedGene);
    if (!g) {
      wrap.innerHTML =
        '<div class="gene-empty"><div class="glyph-lg">\u25C8</div>' +
        "<div><strong>Select a gene</strong> from the panel on the left to see its overview, phenotypes, variants, studies, and publications.</div></div>";
      return;
    }

    const heteroCount = g.variants.filter((v) => v.zygosity === "Heterozygous").length;
    const phasedHetCount = g.variants.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length;
    const phasedPct = heteroCount ? Math.round((phasedHetCount / heteroCount) * 100) : 0;
    const pathCount = g.variants.filter((v) => v.category === "concern").length;
    const protectCount = g.variants.filter((v) => v.category === "protective").length;
    const uncertainCount = g.variants.filter((v) => v.category === "uncertain").length;
    const studyCount = g.variants.reduce((s, v) => s + v.studies.length, 0);

    wrap.innerHTML =
      '<div class="gene-head">' +
        "<h1><span class=\"sym\">" + g.symbol + '</span><span class="full">' + g.name + "</span></h1>" +
        '<button class="btn" id="genome-browser-btn">View in Genome Browser</button>' +
      "</div>" +
      '<div class="gene-coord">' + g.chromosome + " · pLI " + g.pli + " · LOEUF " + g.loeuf + "</div>" +
      '<div class="gene-tabs" id="gene-tabs">' +
        '<button data-tab="overview" class="active">Overview</button>' +
        '<button data-tab="phenotypes">Phenotypes</button>' +
        '<button data-tab="variants">Variants</button>' +
        '<button data-tab="studies">Studies</button>' +
        '<button data-tab="publications">Publications</button>' +
      "</div>" +

      '<div class="gene-pane active" data-pane="overview">' +
        '<div class="kpi-row">' +
          kpi(g.variantsDetected, "Variants detected") +
          kpi(pathCount, "Potential concerns") +
          kpi(protectCount, "Protective") +
          kpi(uncertainCount, "Uncertain") +
          kpi(phasedPct + "%", "Het. variants phased") +
          kpi(g.hpoTermCount, "HPO terms") +
          kpi(g.goTermCount, "GO terms") +
          kpi(g.publications.length, "Publications") +
        "</div>" +
        '<div class="panel-box"><h3>Gene summary</h3><p>' + g.summary + "</p></div>" +
        (g.associatedPathology.length
          ? '<div class="panel-box"><h3>Associated pathology</h3>' +
            g.associatedPathology.map((p) =>
              '<div class="pathology-chip">' + p.name + '<span class="tag">' + p.inheritance + "</span>" +
              (p.omim ? '<span class="tag">OMIM #' + p.omim + "</span>" : "") + "</div>"
            ).join("") +
            "</div>"
          : '<div class="panel-box"><h3>Associated pathology</h3><p style="color:var(--slate);">No single-gene OMIM phenotype is established for ' + g.symbol + "; studied primarily as a quantitative-trait modifier.</p></div>") +
        '<div class="panel-box"><h3>Reference links</h3><div class="ref-links">' +
          '<a class="ref-link" href="' + g.links.ncbiGene + '" target="_blank" rel="noopener">NCBI Gene (' + g.ncbiGeneId + ") \u2197</a>" +
          '<a class="ref-link" href="' + g.links.omim + '" target="_blank" rel="noopener">OMIM *' + g.omimGene + " \u2197</a>" +
          (g.omimPhenotype ? '<a class="ref-link" href="https://omim.org/entry/' + g.omimPhenotype + '" target="_blank" rel="noopener">OMIM phenotype #' + g.omimPhenotype + " \u2197</a>" : "") +
          '<a class="ref-link" href="' + g.links.genecards + '" target="_blank" rel="noopener">GeneCards \u2197</a>' +
          '<a class="ref-link" href="' + g.links.clinvarGene + '" target="_blank" rel="noopener">ClinVar (gene) \u2197</a>' +
        "</div></div>" +
      "</div>" +

      '<div class="gene-pane" data-pane="phenotypes">' +
        (g.hpoTerms.length
          ? '<div class="hpo-card-grid">' + g.hpoTerms.map((t) =>
              '<div class="hpo-card"><div class="id">' + t.id + '</div><div class="label">' + t.label + '</div><div class="evidence">' + t.evidence + "</div></div>"
            ).join("") + "</div>"
          : '<div class="pub-empty">No curated HPO associations recorded for this gene yet.</div>') +
      "</div>" +

      '<div class="gene-pane" data-pane="variants">' + renderGeneVariants(g) + "</div>" +

      '<div class="gene-pane" data-pane="studies">' + renderStudiesTab(g) + "</div>" +

      '<div class="gene-pane" data-pane="publications">' +
        (g.publications.length
          ? '<div class="pub-grid">' + g.publications.map(pubCard).join("") + "</div>"
          : '<div class="pub-empty">No curated publications linked to this gene\u2019s variants yet.</div>') +
      "</div>";

    // Tab switching within gene detail
    document.querySelectorAll("#gene-tabs button").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.querySelectorAll("#gene-tabs button").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        document.querySelectorAll(".gene-pane").forEach((p) => p.classList.remove("active"));
        wrap.querySelector('.gene-pane[data-pane="' + btn.dataset.tab + '"]').classList.add("active");
      });
    });

    byId("genome-browser-btn").addEventListener("click", () => openGenomeModal(g.symbol, g.variants[0]));

    wireVariantsPane(wrap, g);

    // Jump from a Studies card back to the Variants tab for that variant
    wrap.querySelectorAll(".study-card__variant").forEach((el) => {
      el.addEventListener("click", () => {
        const key = g.symbol + ":" + el.dataset.vid;
        state.expandedVariantRows.add(key);
        renderGeneDetail();
        document.querySelector('#gene-tabs button[data-tab="variants"]').click();
      });
    });
  }

  function wireVariantsPane(wrap, g) {
    wrap.querySelectorAll(".vid").forEach((el) => {
      el.addEventListener("click", () => {
        const key = g.symbol + ":" + el.dataset.vid;
        if (state.expandedVariantRows.has(key)) state.expandedVariantRows.delete(key);
        else state.expandedVariantRows.add(key);
        renderGeneDetail();
        document.querySelector('#gene-tabs button[data-tab="variants"]').click();
      });
    });
    wrap.querySelectorAll(".variant-section__head").forEach((el) => {
      el.addEventListener("click", () => {
        const body = el.parentElement.querySelector("tbody");
        const collapsed = body.style.display === "none";
        body.style.display = collapsed ? "" : "none";
        el.querySelector(".caret").textContent = collapsed ? "\u25BE" : "\u25B8";
      });
    });
    wrap.querySelectorAll(".open-genome-modal").forEach((el) => {
      el.addEventListener("click", (e) => {
        e.stopPropagation();
        openGenomeModal(g.symbol, g.variants.find((v) => v.id === el.dataset.vid));
      });
    });
  }

  function kpi(num, label) {
    return '<div class="kpi"><div class="num">' + num + '</div><div class="lbl">' + label + "</div></div>";
  }

  function pubCard(p) {
    return (
      '<div class="pub-card"><h4>' + p.title + "</h4>" +
      '<div class="byline">' + p.authors + " — " + p.journal + ", " + p.year + "</div>" +
      '<div class="finding">' + p.finding + "</div>" +
      '<div class="tags">' + p.tags.map((t) => '<span class="chip">' + t + "</span>").join("") + "</div>" +
      '<a href="https://doi.org/' + p.doi + '" target="_blank" rel="noopener">View on PubMed / DOI \u2197</a>' +
      "</div>"
    );
  }

  function evidenceDots(level) {
    let out = '<span class="evidence-dots" title="Evidence level ' + level + '/3">';
    for (let i = 1; i <= 3; i++) out += '<span class="' + (i <= level ? "filled" : "") + '"></span>';
    return out + "</span>";
  }

  // Aggregates every study across every variant of a gene into its own
  // reviewable list — distinct from the curated Publications tab.
  function renderStudiesTab(g) {
    const rows = [];
    g.variants.forEach((v) => v.studies.forEach((s) => rows.push(Object.assign({ variantId: v.id }, s))));
    if (!rows.length) return '<div class="pub-empty">No individual study findings recorded for this gene\u2019s variants yet.</div>';
    return rows.map((s) =>
      '<div class="study-card">' +
      '<div class="study-card__head"><span class="study-card__variant" data-vid="' + s.variantId + '">' + s.variantId + "</span>" +
      evidenceDots(s.evidenceLevel) + "</div>" +
      '<div class="study-card__condition">' + s.condition + "</div>" +
      '<div class="study-card__finding">' + s.finding + "</div>" +
      '<div class="study-card__foot"><span>' + s.genotypeRelevance + "</span><span>" + (s.source || "") + "</span></div>" +
      "</div>"
    ).join("");
  }

  const CATEGORY_META = {
    concern: { label: "Potential Concerns", cls: "concern" },
    protect: { label: "Protective Associations", cls: "protect" },
    protective: { label: "Protective Associations", cls: "protect" },
    uncertain: { label: "Uncertain Findings", cls: "uncertain" },
    uncategorized: { label: "Not Categorized", cls: "neutral" }
  };

  // -----------------------------------------------------------------
  // Score cell renderers — shared by the gene-level Variants tab and
  // the global Variants view so the two stay in sync.
  // -----------------------------------------------------------------
  function scoreCell(value, max, hotThreshold) {
    if (value === null || value === undefined) return '<span class="na-cell">\u2014</span>';
    const pct = Math.max(0, Math.min(100, (value / max) * 100));
    const hot = hotThreshold !== undefined && value >= hotThreshold;
    return (
      '<div class="score-cell"><div class="score-bar-track"><div class="score-bar-fill' + (hot ? " hot" : "") +
      '" style="width:' + pct + '%"></div></div>' + value + "</div>"
    );
  }

  function readsCell(reads) {
    if (!reads) return '<span class="na-cell">\u2014</span>';
    const pct = Math.round((reads.matching / reads.total) * 100);
    return '<span class="reads-frac">' + reads.matching + "/" + reads.total + '</span> <span class="reads-pct">(' + pct + "%)</span>";
  }

  function scoreColumnsHead() {
    return "<th>CADD</th><th>SpliceAI</th><th>AlphaMissense</th><th>QUAL</th><th>Reads (alt/total)</th>";
  }
  function scoreColumnsCells(v) {
    return (
      "<td>" + scoreCell(v.cadd, 40, 20) + "</td>" +
      "<td>" + scoreCell(v.spliceai, 1, 0.5) + "</td>" +
      "<td>" + (v.alphamissense === null || v.alphamissense === undefined ? '<span class="na-cell">\u2014</span>' : scoreCell(v.alphamissense, 1, 0.564)) + "</td>" +
      "<td>" + (v.qual === null || v.qual === undefined ? '<span class="na-cell">\u2014</span>' : v.qual) + "</td>" +
      "<td>" + readsCell(v.reads) + "</td>"
    );
  }

  function renderGeneVariants(g) {
    const order = ["concern", "protective", "uncertain", "uncategorized"];
    const summary =
      '<div class="cat-summary">' +
      order.map((cat) => {
        const n = g.variants.filter((v) => v.category === cat).length;
        const meta = CATEGORY_META[cat];
        return '<span class="cat-pill cat-pill--' + meta.cls + '"><span class="n">' + n + "</span>" + meta.label + "</span>";
      }).join("") +
      "</div>";

    const sections = order
      .map((cat) => {
        const vs = g.variants.filter((v) => v.category === cat);
        if (!vs.length) return "";
        const meta = CATEGORY_META[cat];
        return (
          '<div class="variant-section">' +
          '<div class="variant-section__head ' + meta.cls + '"><span>' + vs.length + " " + meta.label.toUpperCase() + "</span><span class=\"caret\">\u25BE</span></div>" +
          '<table class="variant-table"><thead><tr>' +
          "<th>Variant ID</th><th>Genotype</th><th>Zygosity</th><th>Phase</th><th>MAF</th><th>Consequence</th><th>ClinVar</th><th>REVEL</th>" +
          scoreColumnsHead() +
          "</tr></thead><tbody>" +
          vs.map((v) => variantRow(g.symbol, v)).join("") +
          "</tbody></table>" +
          "</div>"
        );
      })
      .join("");

    return summary + sections;
  }

  function phaseTag(phase) {
    const cls = (phase || "unknown").toLowerCase().replace(/[^a-z]/g, "");
    const glyph = phase === "Maternal" ? "M" : phase === "Paternal" ? "P" : phase === "N/A" ? "\u2014" : "?";
    return '<span class="phase-tag ' + cls + '"><span class="glyph">[' + glyph + "]</span>" + phase + "</span>";
  }

  function variantRow(geneSymbol, v) {
    const key = geneSymbol + ":" + v.id;
    const expanded = state.expandedVariantRows.has(key);
    let html =
      "<tr>" +
      '<td class="vid" data-vid="' + v.id + '">' + v.id + " " + (expanded ? "\u25BE" : "\u25B8") + "</td>" +
      '<td class="geno">' + v.genotype + "</td>" +
      "<td>" + v.zygosity + "</td>" +
      "<td>" + phaseTag(v.phase) + "</td>" +
      "<td>" + v.maf + "</td>" +
      "<td>" + v.consequence.map((c) => '<span class="chip' + (c === "Missense" ? " chip--missense" : "") + '">' + c + "</span>").join("") + "</td>" +
      "<td>" + v.clinvar + "</td>" +
      "<td>" + (v.revel === null || v.revel === undefined ? '<span class="na-cell">\u2014</span>' : v.revel) + "</td>" +
      scoreColumnsCells(v) +
      "</tr>";

    if (expanded) {
      html +=
        '<tr class="study-row"><td colspan="13">' +
        '<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;gap:10px;flex-wrap:wrap;">' +
        '<div class="ref-links" style="margin-top:0;">' +
        '<a class="ref-link" href="https://www.ncbi.nlm.nih.gov/gene/' + geneBySymbol(geneSymbol).ncbiGeneId + '" target="_blank" rel="noopener">NCBI Gene \u2197</a>' +
        '<a class="ref-link" href="https://omim.org/entry/' + geneBySymbol(geneSymbol).omimGene + '" target="_blank" rel="noopener">OMIM \u2197</a>' +
        "</div>" +
        '<button class="btn open-genome-modal" data-vid="' + v.id + '">View in Genome Browser \u2192 ' + v.coordinate + "</button></div>" +
        (v.studies.length
          ? '<div class="study-line" style="font-weight:700;color:var(--slate);"><div>What studies say</div><div>Condition</div><div>Genotype relevance</div></div>' +
            v.studies.map((s) =>
              '<div class="study-line"><div class="dot-lead">' + s.finding + "</div><div>" + s.condition + '</div><div class="relevance">' + s.genotypeRelevance.replace(/(fully matches|carries the protective|heterozygous carrier|confers deficiency|is heterozygous)/i, "<b>$1</b>") + "</div></div>"
            ).join("")
          : '<div class="pub-empty">Research available; no curated study summary yet.</div>') +
        (v.lastEvaluated ? '<div style="margin-top:8px;font-size:11.5px;color:var(--slate);">Last evaluated: ' + v.lastEvaluated + "</div>" : "") +
        "</td></tr>";
    }
    return html;
  }

  // -----------------------------------------------------------------
  // Genome browser modal (stub / future IGV.js hook)
  // -----------------------------------------------------------------
  function openGenomeModal(symbol, variant) {
    byId("genome-modal-title").textContent = "View in Genome Browser \u2014 " + symbol;
    byId("genome-modal-locus").textContent =
      variant ? variant.coordinate + "  (" + variant.id + ", " + variant.genotype + ")" : "No coordinate available";
    byId("genome-modal").classList.add("open");
  }
  function initGenomeModal() {
    byId("genome-modal-close").addEventListener("click", () => byId("genome-modal").classList.remove("open"));
    byId("genome-modal").addEventListener("click", (e) => {
      if (e.target.id === "genome-modal") byId("genome-modal").classList.remove("open");
    });
  }

  // -----------------------------------------------------------------
  // GENES VIEW
  // -----------------------------------------------------------------
  function computePhasedPct(g) {
    const hetero = g.variants.filter((v) => v.zygosity === "Heterozygous");
    if (!hetero.length) return 0;
    const phased = hetero.filter((v) => v.phase && v.phase !== "Unknown");
    return Math.round((phased.length / hetero.length) * 100);
  }

  let genesSortKey = "symbol", genesSortDir = 1;

  function renderGenesTable() {
    const q = byId("genes-search").value.trim().toLowerCase();
    let rows = GENES.map((g) => ({
      symbol: g.symbol,
      organSystem: g.organSystem,
      variantsDetected: g.variantsDetected,
      pathogenic: g.variants.filter((v) => v.category === "concern").length,
      phasedPct: computePhasedPct(g),
      hpoTermCount: g.hpoTermCount,
      goTermCount: g.goTermCount
    })).filter((r) => !q || r.symbol.toLowerCase().includes(q));

    rows.sort((a, b) => {
      const av = a[genesSortKey], bv = b[genesSortKey];
      if (typeof av === "string") return av.localeCompare(bv) * genesSortDir;
      return (av - bv) * genesSortDir;
    });

    byId("genes-table").querySelector("tbody").innerHTML = rows.map((r) =>
      "<tr data-symbol=\"" + r.symbol + "\">" +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);">' + r.symbol + "</td>" +
      "<td>" + r.organSystem + "</td>" +
      "<td>" + r.variantsDetected + "</td>" +
      "<td>" + (r.pathogenic ? '<span class="badge badge--concern">' + r.pathogenic + "</span>" : "0") + "</td>" +
      '<td><div class="bar-cell"><div class="bar-track"><div class="bar-fill" style="width:' + r.phasedPct + '%"></div></div>' + r.phasedPct + "%</div></td>" +
      "<td>" + r.hpoTermCount + "</td>" +
      "<td>" + r.goTermCount + "</td>" +
      "</tr>"
    ).join("");

    byId("genes-table").querySelectorAll("tbody tr").forEach((tr) => {
      tr.addEventListener("click", () => {
        document.querySelector('.tabs button[data-view="ontology"]').click();
        selectGene(tr.dataset.symbol);
      });
    });
  }

  function initGenesView() {
    byId("genes-search").addEventListener("input", renderGenesTable);
    byId("genes-table").querySelectorAll("th[data-sort]").forEach((th) => {
      th.addEventListener("click", () => {
        const key = th.dataset.sort;
        genesSortDir = genesSortKey === key ? -genesSortDir : 1;
        genesSortKey = key;
        renderGenesTable();
      });
    });
    renderGenesTable();
  }

  // -----------------------------------------------------------------
  // VARIANTS VIEW
  // -----------------------------------------------------------------
  function allVariantsFlat() {
    const rows = [];
    GENES.forEach((g) => g.variants.forEach((v) => rows.push(Object.assign({ gene: g.symbol }, v))));
    return rows;
  }

  function renderVariantsTable() {
    const geneF = byId("filter-gene").value;
    const clinvarF = byId("filter-clinvar").value;
    const phaseF = byId("filter-phase").value;
    const zygF = byId("filter-zygosity").value;

    const rows = allVariantsFlat().filter((v) =>
      (!geneF || v.gene === geneF) &&
      (!clinvarF || v.clinvar === clinvarF) &&
      (!phaseF || v.phase === phaseF) &&
      (!zygF || v.zygosity === zygF)
    );

    const badgeCls = { concern: "badge--concern", protective: "badge--protect", uncertain: "badge--uncertain", uncategorized: "badge--neutral" };

    byId("variants-table").querySelector("tbody").innerHTML = rows.map((v) =>
      "<tr data-gene=\"" + v.gene + "\" data-vid=\"" + v.id + "\">" +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);">' + v.id + "</td>" +
      '<td class="mono">' + v.gene + "</td>" +
      "<td><span class=\"badge " + (badgeCls[v.category] || "badge--neutral") + "\">" + v.clinvar + "</span></td>" +
      "<td>" + (v.revel === null || v.revel === undefined ? '<span class="na-cell">\u2014</span>' : v.revel) + "</td>" +
      scoreColumnsCells(v) +
      '<td class="mono">' + v.coordinate + "</td>" +
      "<td>" + v.zygosity + "</td>" +
      "<td>" + phaseTag(v.phase) + "</td>" +
      "<td>" + v.maf + "</td>" +
      "<td>" + v.consequence.map((c) => '<span class="chip">' + c + "</span>").join("") + "</td>" +
      "</tr>"
    ).join("");

    byId("variants-table").querySelectorAll("tbody tr").forEach((tr) => {
      tr.addEventListener("click", () => {
        const vid = tr.dataset.vid;
        document.querySelector('.tabs button[data-view="ontology"]').click();
        state.expandedVariantRows.add(tr.dataset.gene + ":" + vid);
        selectGene(tr.dataset.gene);
        document.querySelector('#gene-tabs button[data-tab="variants"]').click();
      });
    });
  }

  function initVariantsView() {
    const geneSel = byId("filter-gene");
    GENES.forEach((g) => {
      const opt = document.createElement("option");
      opt.value = g.symbol;
      opt.textContent = g.symbol;
      geneSel.appendChild(opt);
    });
    ["filter-gene", "filter-clinvar", "filter-phase", "filter-zygosity"].forEach((id) =>
      byId(id).addEventListener("change", renderVariantsTable)
    );
    renderVariantsTable();
  }

  // -----------------------------------------------------------------
  // ANALYSIS VIEW
  // -----------------------------------------------------------------
  function renderAnalysisView() {
    const totalVariants = GENES.reduce((s, g) => s + g.variantsDetected, 0);
    const totalPath = GENES.reduce((s, g) => s + g.variants.filter((v) => v.category === "concern").length, 0);
    const totalHetero = GENES.reduce((s, g) => s + g.variants.filter((v) => v.zygosity === "Heterozygous").length, 0);
    const totalPhased = GENES.reduce((s, g) => s + g.variants.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length, 0);
    const phasedPct = totalHetero ? Math.round((totalPhased / totalHetero) * 100) : 0;

    byId("analysis-kpis").innerHTML =
      kpi(totalVariants.toLocaleString(), "Total variants (curated genes)") +
      kpi(totalPath, "Pathogenic / LP calls") +
      kpi(phasedPct + "%", "Heterozygous calls phased") +
      kpi(GENES.length, "Genes in panel");

    byId("prs-grid").innerHTML = PRS.map((p) =>
      '<div class="prs-card"><div class="prs-card__head"><span class="trait">' + p.trait + '</span><span class="badge badge--' +
      (p.category === "HIGH" ? "concern" : p.category === "PROTECTIVE" ? "protect" : "uncertain") + '">' + p.category + "</span></div>" +
      '<div class="prs-track"><div class="prs-fill ' + p.category + '" style="width:' + p.percentile + '%"></div></div>' +
      '<div class="prs-foot"><span>' + p.percentile + "th percentile · " + p.organSystem + "</span><a href=\"#\">" + p.pgsId + "</a></div>" +
      "</div>"
    ).join("");

    byId("pgx-body").innerHTML = PGX.map((p) =>
      "<tr><td class=\"mono\" style=\"font-weight:700;\">" + p.gene + '</td><td class="mono">' + p.diplotype + "</td><td>" + p.phenotype +
      "</td><td>" + p.drug + '</td><td><span class="action-tier ' + p.actionTier + '">' + p.actionTier + "</span></td><td>" + p.recommendation + "</td></tr>"
    ).join("");
  }

  // -----------------------------------------------------------------
  // REPORTS VIEW
  // -----------------------------------------------------------------
  function renderReportsView() {
    byId("report-meta").textContent = REPORT.sampleLabel + " · Generated " + REPORT.generated;
    byId("report-narrative").innerHTML = "<p>" + REPORT.narrative + "</p>";
    byId("report-gene-body").innerHTML = REPORT.geneBreakdown.map((r) =>
      "<tr><td class=\"mono\" style=\"font-weight:700;color:var(--teal-dark);\">" + r.symbol + "</td><td>" + r.variantsDetected +
      "</td><td>" + (r.pathogenicOrLP ? '<span class="badge badge--concern">' + r.pathogenicOrLP + "</span>" : "0") +
      "</td><td>" + (r.protective ? '<span class="badge badge--protect">' + r.protective + "</span>" : "0") +
      "</td><td>" + (r.uncertain ? '<span class="badge badge--uncertain">' + r.uncertain + "</span>" : "0") + "</td></tr>"
    ).join("");

    byId("download-json-btn").addEventListener("click", () => {
      const payload = { jobMeta: JOB_META, report: REPORT, genes: GENES, prs: PRS, pgx: PGX };
      const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "genomic-report.json";
      a.click();
      URL.revokeObjectURL(url);
    });
  }

  // -----------------------------------------------------------------
  // Decorative DNA watermark (colorful double helix, faint, generated
  // once at load — see .dna-watermark in css/style.css for how it's
  // kept out from under actual data)
  // -----------------------------------------------------------------
  function buildDnaWatermark() {
    const g = byId("dna-rungs");
    if (!g) return;
    const svgns = "http://www.w3.org/2000/svg";
    const height = 900, amplitude = 70, centerX = 200, period = 220, rungGap = 30;
    const baseColor = { A: "#1e7b4d", T: "#b3261e", G: "#2a5c99", C: "#b5722a" };
    const bases = ["A", "T", "G", "C"];

    let dA = "", dB = "";
    const steps = 90;
    for (let i = 0; i <= steps; i++) {
      const y = (height / steps) * i;
      const angle = (y / period) * 2 * Math.PI;
      const xA = centerX + amplitude * Math.sin(angle);
      const xB = centerX - amplitude * Math.sin(angle);
      dA += (i === 0 ? "M" : "L") + xA.toFixed(1) + "," + y.toFixed(1) + " ";
      dB += (i === 0 ? "M" : "L") + xB.toFixed(1) + "," + y.toFixed(1) + " ";
    }
    const pathA = document.createElementNS(svgns, "path");
    pathA.setAttribute("d", dA); pathA.setAttribute("class", "strand"); pathA.setAttribute("stroke", "#1c5b5e");
    const pathB = document.createElementNS(svgns, "path");
    pathB.setAttribute("d", dB); pathB.setAttribute("class", "strand"); pathB.setAttribute("stroke", "#2a5c99");
    g.appendChild(pathA);
    g.appendChild(pathB);

    let idx = 0;
    for (let y = 10; y < height; y += rungGap) {
      const angle = (y / period) * 2 * Math.PI;
      const xA = centerX + amplitude * Math.sin(angle);
      const xB = centerX - amplitude * Math.sin(angle);
      const base = bases[idx % 4]; idx++;
      const line = document.createElementNS(svgns, "line");
      line.setAttribute("x1", xA.toFixed(1)); line.setAttribute("y1", y);
      line.setAttribute("x2", xB.toFixed(1)); line.setAttribute("y2", y);
      line.setAttribute("class", "rung"); line.setAttribute("stroke", baseColor[base]);
      g.appendChild(line);

      const label = document.createElementNS(svgns, "text");
      label.setAttribute("x", (Math.max(xA, xB) + 6).toFixed(1));
      label.setAttribute("y", (y + 4).toFixed(1));
      label.setAttribute("fill", baseColor[base]);
      label.textContent = base;
      g.appendChild(label);
    }
  }

  // -----------------------------------------------------------------
  // Init
  // -----------------------------------------------------------------
  document.addEventListener("DOMContentLoaded", () => {
    buildDnaWatermark();
    initTabs();
    initOntologySwitch();
    initGenomeModal();
    renderLeftPanel();
    renderGeneDetail();
    initGenesView();
    initVariantsView();
    renderAnalysisView();
    renderReportsView();
  });
})();
