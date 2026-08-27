/**
 * Genomic Ontology Explorer — app logic (v3)
 * Full Multi-Level Hierarchy & Inspection (Level 1, 2, 3, 4)
 * Real OpenCRAVAT dataset integration with Genome Browser lollipop visualizer.
 */

(function () {
  "use strict";

  // -----------------------------------------------------------------
  // Shared state
  // -----------------------------------------------------------------
  const state = {
    tree: "hpo",
    layout: "tree",               // default to Tree for rich 4-level navigation
    scopeFindingsOnly: false,
    treeSearch: "",
    selectedTarget: null,         // { type: 'gene'|'category', symbol?, id?, label?, level?, genes:[] }
    expandedNodes: new Set(),      // ontology group/term ids currently expanded
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
        const targetView = byId("view-" + btn.dataset.view);
        if (targetView) targetView.classList.add("active");

        if (btn.dataset.view === "genes") renderGenesTable();
        else if (btn.dataset.view === "variants") renderVariantsTable();
        else if (btn.dataset.view === "analysis") renderAnalysisView();
        else if (btn.dataset.view === "reports") renderReportsView();
      });
    });

    if (window.JOB_META) {
      byId("job-meta").textContent =
        JOB_META.sample + " · OpenCRAVAT " + JOB_META.opencravatVersion + " · " +
        (JOB_META.uniqueVariants || 0).toLocaleString() + " variants";
    }
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

    const scopeCheck = byId("scope-findings-only");
    if (scopeCheck) {
      scopeCheck.addEventListener("change", (e) => {
        state.scopeFindingsOnly = e.target.checked;
        renderLeftPanel();
      });
    }

    const searchInput = byId("tree-search");
    if (searchInput) {
      searchInput.addEventListener("input", (e) => {
        state.treeSearch = e.target.value.trim().toLowerCase();
        renderLeftPanel();
      });
    }
  }

  function matchesSearch(text) {
    if (!state.treeSearch) return true;
    return (text || "").toLowerCase().includes(state.treeSearch);
  }

  function geneHasFindings(symbol) {
    const g = geneBySymbol(symbol);
    if (!g) return false;
    return g.variants.some((v) => v.category === "concern" || v.category === "protective" || v.category === "uncertain");
  }

  // -----------------------------------------------------------------
  // Filtered Ontology Data Model (Unified 4-Level)
  // -----------------------------------------------------------------
  function filteredOntology() {
    const ont = ONTOLOGIES[state.tree] || ONTOLOGIES["hpo"];
    if (!ont || !ont.groups) return [];
    const out = [];

    ont.groups.forEach((group) => {
      let groupGenes = Array.from(new Set(group.genes || []));
      if (state.scopeFindingsOnly) groupGenes = groupGenes.filter(geneHasFindings);
      if (state.scopeFindingsOnly && !groupGenes.length) return;

      const subcategories = (group.terms || []).map((subcat) => {
        let subGenes = Array.from(new Set(subcat.genes || []));
        if (state.scopeFindingsOnly) subGenes = subGenes.filter(geneHasFindings);
        if (state.scopeFindingsOnly && !subGenes.length) return null;

        const terms = (subcat.terms || []).map((term) => {
          let termGenes = Array.from(new Set(term.genes || []));
          if (state.scopeFindingsOnly) termGenes = termGenes.filter(geneHasFindings);
          if (state.treeSearch) termGenes = termGenes.filter((s) => matchesSearch(s) || matchesSearch(term.label) || matchesSearch(term.id));
          if (state.treeSearch && !termGenes.length && !matchesSearch(term.label) && !matchesSearch(term.id)) return null;
          return { term, genes: termGenes };
        }).filter(Boolean);

        const subcatHit = matchesSearch(subcat.label) || matchesSearch(subcat.id) || subGenes.some(matchesSearch) || terms.length > 0;
        if (state.treeSearch && !subcatHit) return null;

        return { subcategory: subcat, genes: subGenes, terms };
      }).filter(Boolean);

      const groupHit = matchesSearch(group.label) || matchesSearch(group.id) || groupGenes.some(matchesSearch) || subcategories.length > 0;
      if (state.treeSearch && !groupHit) return;

      out.push({ group, genes: groupGenes, subcategories });
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
  // Tree / List Rendering with Multi-Level Selection
  // -----------------------------------------------------------------
  function renderTreeOrList() {
    const data = filteredOntology();
    const scroll = byId("tree-scroll");
    scroll.className = "tree-scroll" + (state.layout === "list" ? " tree-list-mode" : "");
    scroll.innerHTML = "";

    data.forEach(({ group, genes: groupGenes, subcategories }) => {
      const totalGenes = groupGenes.length;
      const groupEl = document.createElement("div");
      groupEl.className = "tree-group";

      const isOpen = state.expandedNodes.has(group.id) || !!state.treeSearch;
      const isSelected = state.selectedTarget && state.selectedTarget.id === group.id;

      // LEVEL 1: System / Root Group
      const row = document.createElement("div");
      row.className = "tree-row tree-row--organ" + (isOpen ? " expanded" : "") + (isSelected ? " selected" : "");
      row.innerHTML =
        '<span class="caret">' + (isOpen ? "\u25BE" : "\u25B8") + '</span><span class="node-dot"></span>' +
        "<span>" + group.label + "</span>" +
        '<span class="count-chip">' + totalGenes + " gene" + (totalGenes === 1 ? "" : "s") + "</span>";

      row.addEventListener("click", (e) => {
        // Toggle expansion
        if (state.expandedNodes.has(group.id)) state.expandedNodes.delete(group.id);
        else state.expandedNodes.add(group.id);

        // Select Level 1 Category
        selectCategory({
          type: "category",
          level: 1,
          id: group.id,
          label: group.label,
          levelName: "Organ System / Root Category",
          genes: groupGenes,
          subcategories
        });
      });
      groupEl.appendChild(row);

      const childWrap = document.createElement("div");
      childWrap.className = "tree-children" + (isOpen ? " open" : "");

      if (state.layout === "tree" && subcategories.length) {
        subcategories.forEach(({ subcategory: subcat, genes: subGenes, terms }) => {
          const subOpen = state.expandedNodes.has(subcat.id) || !!state.treeSearch;
          const subSelected = state.selectedTarget && state.selectedTarget.id === subcat.id;

          // LEVEL 2: Subcategory
          const subRow = document.createElement("div");
          subRow.className = "tree-row tree-row--subcat" + (subOpen ? " expanded" : "") + (subSelected ? " selected" : "");
          subRow.innerHTML =
            '<span class="caret">' + (subOpen ? "\u25BE" : "\u25B8") + '</span><span class="node-dot" style="background:var(--teal);"></span>' +
            "<span>" + subcat.label + "</span>" +
            '<span class="count-chip">' + subGenes.length + "</span>";

          subRow.addEventListener("click", (e) => {
            e.stopPropagation();
            if (state.expandedNodes.has(subcat.id)) state.expandedNodes.delete(subcat.id);
            else state.expandedNodes.add(subcat.id);

            selectCategory({
              type: "category",
              level: 2,
              id: subcat.id,
              label: subcat.label,
              levelName: "Subcategory / Clinical Partition",
              genes: subGenes,
              parentLabel: group.label,
              terms
            });
          });
          childWrap.appendChild(subRow);

          const subChildWrap = document.createElement("div");
          subChildWrap.className = "tree-children" + (subOpen ? " open" : "");

          if (terms && terms.length) {
            terms.forEach(({ term, genes: termGenes }) => {
              const termOpen = state.expandedNodes.has(term.id) || !!state.treeSearch;
              const termSelected = state.selectedTarget && state.selectedTarget.id === term.id;

              // LEVEL 3: Phenotype / GO Term
              const termRow = document.createElement("div");
              termRow.className = "tree-row tree-row--term" + (termOpen ? " expanded" : "") + (termSelected ? " selected" : "");
              termRow.innerHTML =
                '<span class="caret">' + (termOpen ? "\u25BE" : "\u25B8") + '</span><span class="node-dot"></span>' +
                "<span>" + term.label + '</span><span class="count-chip mono">' + term.id + "</span>";

              termRow.addEventListener("click", (e) => {
                e.stopPropagation();
                if (state.expandedNodes.has(term.id)) state.expandedNodes.delete(term.id);
                else state.expandedNodes.add(term.id);

                selectCategory({
                  type: "category",
                  level: 3,
                  id: term.id,
                  label: term.label,
                  levelName: "Phenotype / Functional Term (" + term.id + ")",
                  genes: termGenes,
                  parentLabel: subcat.label
                });
              });
              subChildWrap.appendChild(termRow);

              // LEVEL 4: Genes under this Term
              const geneWrap = document.createElement("div");
              geneWrap.className = "tree-children" + (termOpen ? " open" : "");
              termGenes.forEach((sym) => geneWrap.appendChild(makeGeneRow(sym)));
              subChildWrap.appendChild(geneWrap);
            });
          } else {
            // Direct genes under subcategory
            subGenes.forEach((sym) => subChildWrap.appendChild(makeGeneRow(sym)));
          }

          childWrap.appendChild(subChildWrap);
        });
      } else {
        // Flat List mode: System -> All Genes directly
        groupGenes.forEach((sym) => childWrap.appendChild(makeGeneRow(sym)));
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
    const isSelected = state.selectedTarget && state.selectedTarget.type === "gene" && state.selectedTarget.symbol === symbol;
    row.className = "tree-row tree-row--gene" + (isSelected ? " selected" : "");
    const flagged = g && g.variants.some((v) => v.category === "concern");
    row.innerHTML =
      '<span class="node-dot"></span><span>' + symbol + "</span>" +
      (flagged ? '<span class="count-chip" style="color:var(--concern);border-color:var(--concern-bg);font-weight:700;">flag</span>' : "");
    row.addEventListener("click", (e) => {
      e.stopPropagation();
      selectGene(symbol);
    });
    return row;
  }

  function selectGene(symbol) {
    state.selectedTarget = { type: "gene", symbol };
    if (state.layout === "graph") renderGraph(); else renderTreeOrList();
    renderInspectionPanel();
  }

  function selectCategory(catObj) {
    state.selectedTarget = catObj;
    if (state.layout === "graph") renderGraph(); else renderTreeOrList();
    renderInspectionPanel();
  }

  // -----------------------------------------------------------------
  // Graph layout mode — Interactive Multi-Level Node/Edge Diagram
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
    const ROW_H = 32, GROUP_X = 26, SUBCAT_X = 170, GENE_X = 340, PAD_TOP = 20, WIDTH = 680;

    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    svg.setAttribute("class", "graph-svg");
    svg.setAttribute("width", WIDTH);

    const geneChips = []; // { y, genes:[] }
    let y = PAD_TOP;
    const edges = [];
    const nodes = [];

    data.forEach(({ group, genes: groupGenes, subcategories }) => {
      const groupY = y;
      nodes.push({ type: "group", id: group.id, x: GROUP_X, y: groupY, label: group.label, genes: groupGenes, level: 1 });
      y += ROW_H;

      subcategories.forEach(({ subcategory: subcat, genes: subGenes }) => {
        const subcatY = y;
        nodes.push({ type: "subcat", id: subcat.id, x: SUBCAT_X, y: subcatY, label: subcat.label, genes: subGenes, level: 2 });
        edges.push({ x1: GROUP_X + 6, y1: groupY, x2: SUBCAT_X - 6, y2: subcatY });
        geneChips.push({ y: subcatY, genes: subGenes });
        y += ROW_H;
      });
      y += 8; // gap between groups
    });

    const totalHeight = y + 20;
    svg.setAttribute("height", totalHeight);
    svg.setAttribute("viewBox", "0 0 " + WIDTH + " " + totalHeight);

    // Draw Edges
    edges.forEach((e) => {
      const line = document.createElementNS(svgns, "line");
      line.setAttribute("x1", e.x1); line.setAttribute("y1", e.y1);
      line.setAttribute("x2", e.x2); line.setAttribute("y2", e.y2);
      line.setAttribute("class", "graph-edge");
      line.setAttribute("stroke", "var(--line)");
      line.setAttribute("stroke-width", "1.5");
      svg.appendChild(line);
    });

    // Draw Nodes
    nodes.forEach((n) => {
      const isSel = state.selectedTarget && state.selectedTarget.id === n.id;
      const gNode = document.createElementNS(svgns, "g");
      gNode.style.cursor = "pointer";

      const circle = document.createElementNS(svgns, "circle");
      circle.setAttribute("cx", n.x); circle.setAttribute("cy", n.y);
      circle.setAttribute("r", n.type === "group" ? 6 : 4.5);
      circle.setAttribute("fill", isSel ? "var(--concern)" : (n.type === "group" ? "var(--teal-dark)" : "var(--teal)"));
      gNode.appendChild(circle);

      const text = document.createElementNS(svgns, "text");
      text.setAttribute("x", n.x + 10);
      text.setAttribute("y", n.y + 4);
      text.setAttribute("class", "graph-node-label " + n.type);
      text.setAttribute("font-size", n.type === "group" ? "12px" : "11px");
      text.setAttribute("font-weight", n.type === "group" ? "700" : "600");
      text.setAttribute("fill", isSel ? "var(--teal-dark)" : "var(--ink)");
      text.textContent = n.label.length > 24 ? n.label.substring(0, 22) + "…" : n.label;
      gNode.appendChild(text);

      gNode.addEventListener("click", () => {
        selectCategory({
          type: "category",
          level: n.level,
          id: n.id,
          label: n.label,
          levelName: n.type === "group" ? "Organ System / Root Category" : "Subcategory",
          genes: n.genes
        });
      });

      svg.appendChild(gNode);
    });

    host.appendChild(svg);

    // Gene chips layer
    const geneLayer = document.createElement("div");
    geneLayer.className = "graph-genes-layer";
    geneLayer.style.width = WIDTH + "px";
    geneLayer.style.height = totalHeight + "px";
    geneChips.forEach((row) => {
      row.genes.slice(0, 4).forEach((sym, i) => {
        const chip = document.createElement("div");
        const isGeneSel = state.selectedTarget && state.selectedTarget.type === "gene" && state.selectedTarget.symbol === sym;
        chip.className = "graph-gene-chip" + (isGeneSel ? " selected" : "");
        chip.textContent = sym;
        chip.style.left = (GENE_X + i * 68) + "px";
        chip.style.top = (row.y - 10) + "px";
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
  // Inspection Router: Gene View OR Category Rollup View
  // -----------------------------------------------------------------
  function renderInspectionPanel() {
    const wrap = byId("gene-detail");
    if (!state.selectedTarget) {
      wrap.innerHTML =
        '<div class="gene-empty"><div class="glyph-lg">\u25C8</div>' +
        "<div><strong>Select any category or gene</strong> from the tree on the left to see its rolled-up overview, phenotypes, variants, studies, and publications.</div></div>";
      return;
    }

    if (state.selectedTarget.type === "gene") {
      renderSingleGeneView(wrap, state.selectedTarget.symbol);
    } else {
      renderCategoryRollupView(wrap, state.selectedTarget);
    }
  }

  // -----------------------------------------------------------------
  // 1. Single Gene Detail View
  // -----------------------------------------------------------------
  function renderSingleGeneView(wrap, symbol) {
    const g = geneBySymbol(symbol);
    if (!g) return;

    const heteroCount = g.variants.filter((v) => v.zygosity === "Heterozygous").length;
    const phasedHetCount = g.variants.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length;
    const phasedPct = heteroCount ? Math.round((phasedHetCount / heteroCount) * 100) : 0;
    const pathCount = g.variants.filter((v) => v.category === "concern").length;
    const protectCount = g.variants.filter((v) => v.category === "protective").length;
    const uncertainCount = g.variants.filter((v) => v.category === "uncertain").length;

    wrap.innerHTML =
      '<div class="gene-head">' +
        '<div>' +
          '<h1><span class="sym">' + g.symbol + '</span> <span class="full">' + g.name + '</span></h1>' +
          '<div class="gene-coord mono" style="font-size:11.5px;color:var(--slate);margin-top:2px;">' + g.chromosome + " · pLI " + g.pli + " · LOEUF " + g.loeuf + " · " + g.organSystem + '</div>' +
        '</div>' +
        '<button class="btn btn--primary" id="genome-browser-btn">View in Genome Browser</button>' +
      '</div>' +
      '<div class="gene-tabs" id="gene-tabs">' +
        '<button data-tab="overview" class="active">Overview</button>' +
        '<button data-tab="phenotypes">Phenotypes (' + g.hpoTermCount + ')</button>' +
        '<button data-tab="variants">Variants (' + g.variants.length + ')</button>' +
        '<button data-tab="studies">Studies</button>' +
        '<button data-tab="publications">Publications (' + g.publications.length + ')</button>' +
      '</div>' +

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
        '</div>' +
        '<div class="panel-box"><h3>Gene summary</h3><p>' + g.summary + '</p></div>' +
        (g.associatedPathology && g.associatedPathology.length
          ? '<div class="panel-box"><h3>Associated pathology</h3>' +
            g.associatedPathology.map((p) =>
              '<div class="pathology-chip">' + p.name + '<span class="tag">' + p.inheritance + '</span>' +
              (p.omim ? '<span class="tag">OMIM #' + p.omim + '</span>' : "") + '</div>'
            ).join("") + '</div>'
          : '<div class="panel-box"><h3>Associated pathology</h3><p style="color:var(--slate);">No single-gene OMIM phenotype established; evaluated in quantitative/complex trait networks.</p></div>') +
        '<div class="panel-box"><h3>Reference links</h3><div class="ref-links">' +
          '<a class="ref-link" href="' + g.links.ncbiGene + '" target="_blank" rel="noopener">NCBI Gene (' + g.ncbiGeneId + ') ↗</a>' +
          '<a class="ref-link" href="' + g.links.omim + '" target="_blank" rel="noopener">OMIM *' + g.omimGene + ' ↗</a>' +
          (g.omimPhenotype ? '<a class="ref-link" href="https://omim.org/entry/' + g.omimPhenotype + '" target="_blank" rel="noopener">OMIM Phenotype #' + g.omimPhenotype + ' ↗</a>' : "") +
          '<a class="ref-link" href="' + g.links.genecards + '" target="_blank" rel="noopener">GeneCards ↗</a>' +
          '<a class="ref-link" href="' + g.links.clinvarGene + '" target="_blank" rel="noopener">ClinVar (gene) ↗</a>' +
        '</div></div>' +
      '</div>' +

      '<div class="gene-pane" data-pane="phenotypes">' +
        (g.hpoTerms && g.hpoTerms.length
          ? '<div class="hpo-card-grid">' + g.hpoTerms.map((t) =>
              '<div class="hpo-card"><div class="id">' + t.id + '</div><div class="label">' + t.label + '</div><div class="evidence">' + t.evidence + '</div></div>'
            ).join("") + '</div>'
          : '<div class="pub-empty">No curated HPO associations recorded for this gene.</div>') +
      '</div>' +

      '<div class="gene-pane" data-pane="variants">' + renderVariantSectionList(g.variants, g.symbol) + '</div>' +

      '<div class="gene-pane" data-pane="studies">' + renderStudiesTabForVariants(g.variants) + '</div>' +

      '<div class="gene-pane" data-pane="publications">' +
        (g.publications && g.publications.length
          ? '<div class="pub-grid">' + g.publications.map(pubCard).join("") + '</div>'
          : '<div class="pub-empty">No curated publications directly linked to this gene’s variants yet.</div>') +
      '</div>';

    wireTabSwitching(wrap);
    byId("genome-browser-btn").addEventListener("click", () => openGenomeModal(g.symbol, g.variants));
    wireVariantsPane(wrap);
  }

  // -----------------------------------------------------------------
  // 2. Category / Domain Rollup View (Multi-Level Aggregation)
  // -----------------------------------------------------------------
  function renderCategoryRollupView(wrap, cat) {
    const geneObjects = (cat.genes || []).map(geneBySymbol).filter(Boolean);
    const allVariants = geneObjects.flatMap((g) => g.variants);

    const totalVariants = allVariants.length;
    const pathCount = allVariants.filter((v) => v.category === "concern").length;
    const protectCount = allVariants.filter((v) => v.category === "protective").length;
    const uncertainCount = allVariants.filter((v) => v.category === "uncertain").length;
    const heteroCount = allVariants.filter((v) => v.zygosity === "Heterozygous").length;
    const phasedHetCount = allVariants.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length;
    const phasedPct = heteroCount ? Math.round((phasedHetCount / heteroCount) * 100) : 0;

    // Collect all HPO terms across genes in domain
    const hpoTermsMap = {};
    geneObjects.forEach((g) => {
      (g.hpoTerms || []).forEach((t) => {
        if (!hpoTermsMap[t.id]) hpoTermsMap[t.id] = { id: t.id, label: t.label, genes: [] };
        hpoTermsMap[t.id].genes.push(g.symbol);
      });
    });
    const domainHpos = Object.values(hpoTermsMap);

    wrap.innerHTML =
      '<div class="category-head">' +
        '<div>' +
          '<h1><span>' + cat.label + '</span> <span class="category-level-badge">' + (cat.levelName || "Category") + '</span></h1>' +
          '<div class="category-meta mono">' + (cat.parentLabel ? cat.parentLabel + " \u2192 " : "") + cat.label + " · " + geneObjects.length + " Genes · " + totalVariants + " Actionable Variants</div>" +
        '</div>' +
      '</div>' +

      '<div class="gene-tabs" id="gene-tabs">' +
        '<button data-tab="overview" class="active">Overview</button>' +
        '<button data-tab="phenotypes">Phenotypes (' + domainHpos.length + ')</button>' +
        '<button data-tab="variants">Variants (' + totalVariants + ')</button>' +
        '<button data-tab="studies">Studies</button>' +
        '<button data-tab="publications">Publications</button>' +
      '</div>' +

      '<div class="gene-pane active" data-pane="overview">' +
        '<div class="kpi-row">' +
          kpi(geneObjects.length, "Genes in domain") +
          kpi(pathCount, "Potential concerns") +
          kpi(protectCount, "Protective") +
          kpi(uncertainCount, "Uncertain") +
          kpi(phasedPct + "%", "Het. variants phased") +
          kpi(domainHpos.length, "HPO terms") +
          kpi(totalVariants, "Total variants") +
        '</div>' +

        '<div class="panel-box">' +
          '<h3>Clinical & Biological Overview</h3>' +
          '<p>Roll-up summary of all findings classified under <strong>' + cat.label + '</strong>. This branch comprises <strong>' + geneObjects.length + ' genes</strong> and <strong>' + totalVariants + ' actionable variants</strong> identified in the patient’s phased sequencing data.</p>' +
        '</div>' +

        '<div class="panel-box">' +
          '<h3>Genes in this Branch (' + geneObjects.length + ')</h3>' +
          '<div class="gene-chip-roster">' +
            geneObjects.map((g) => {
              const hasConcern = g.variants.some((v) => v.category === "concern");
              return '<div class="gene-chip-item' + (hasConcern ? " has-concern" : "") + '" data-gene="' + g.symbol + '">' +
                '<span>' + g.symbol + '</span>' +
                '<span class="mono" style="font-size:10px;opacity:0.8;">(' + g.variants.length + 'v)</span>' +
                (hasConcern ? '<span style="font-size:9px;font-weight:700;">●</span>' : '') +
              '</div>';
            }).join("") +
          '</div>' +
        '</div>' +
      '</div>' +

      '<div class="gene-pane" data-pane="phenotypes">' +
        (domainHpos.length
          ? '<div class="hpo-card-grid">' + domainHpos.slice(0, 40).map((t) =>
              '<div class="hpo-card">' +
                '<div class="id">' + t.id + '</div>' +
                '<div class="label">' + t.label + '</div>' +
                '<div class="evidence" style="margin-top:6px;font-size:11px;color:var(--teal-dark);font-weight:600;">Genes: ' + t.genes.join(", ") + '</div>' +
              '</div>'
            ).join("") + '</div>'
          : '<div class="pub-empty">No curated HPO terms assigned to genes in this category.</div>') +
      '</div>' +

      '<div class="gene-pane" data-pane="variants">' +
        renderVariantSectionList(allVariants, cat.label) +
      '</div>' +

      '<div class="gene-pane" data-pane="studies">' +
        renderStudiesTabForVariants(allVariants) +
      '</div>' +

      '<div class="gene-pane" data-pane="publications">' +
        '<div class="pub-grid">' +
          geneObjects.flatMap((g) => g.publications || []).map(pubCard).join("") +
        '</div>' +
        (geneObjects.every((g) => !g.publications || !g.publications.length) ? '<div class="pub-empty">No curated publications for this category.</div>' : "") +
      '</div>';

    wireTabSwitching(wrap);

    // Clicking any gene chip inside the rollup view drills down to that single gene!
    wrap.querySelectorAll(".gene-chip-item").forEach((chip) => {
      chip.addEventListener("click", () => {
        selectGene(chip.dataset.gene);
      });
    });

    wireVariantsPane(wrap);
  }

  function wireTabSwitching(wrap) {
    wrap.querySelectorAll("#gene-tabs button").forEach((btn) => {
      btn.addEventListener("click", () => {
        wrap.querySelectorAll("#gene-tabs button").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        wrap.querySelectorAll(".gene-pane").forEach((p) => p.classList.remove("active"));
        const targetPane = wrap.querySelector('.gene-pane[data-pane="' + btn.dataset.tab + '"]');
        if (targetPane) targetPane.classList.add("active");
      });
    });
  }

  // -----------------------------------------------------------------
  // Variants Section Generator (shared across Gene and Category views)
  // -----------------------------------------------------------------
  const CATEGORY_META = {
    concern: { label: "Potential concerns", cls: "concern" },
    protective: { label: "Protective associations", cls: "protect" },
    uncertain: { label: "Uncertain significance", cls: "uncertain" },
    uncategorized: { label: "Not categorized / Research", cls: "neutral" }
  };

  function renderVariantSectionList(variants, contextSymbol) {
    if (!variants || !variants.length) return '<div class="pub-empty">No variants recorded under this selection.</div>';

    const grouped = { concern: [], protective: [], uncertain: [], uncategorized: [] };
    variants.forEach((v) => {
      const cat = grouped[v.category] ? v.category : "uncategorized";
      grouped[cat].push(v);
    });

    const summary =
      '<div class="kpi-row" style="margin-bottom:14px;">' +
      kpi(variants.length, "Total variants") +
      kpi(grouped.concern.length, "Concerns") +
      kpi(grouped.protective.length, "Protective") +
      kpi(grouped.uncertain.length, "Uncertain") +
      "</div>";

    const sections = ["concern", "protective", "uncertain", "uncategorized"]
      .map((cat) => {
        const vs = grouped[cat];
        if (!vs.length) return "";
        const meta = CATEGORY_META[cat];
        return (
          '<div class="variant-section">' +
          '<div class="variant-section__head ' + meta.cls + '"><span>' + vs.length + " " + meta.label.toUpperCase() + '</span><span class="caret">\u25BE</span></div>' +
          '<div class="data-table-wrap"><table class="variant-table"><thead><tr>' +
          "<th>Variant ID</th><th>Gene</th><th>Genotype</th><th>Zygosity</th><th>Phase</th><th>MAF</th><th>Consequence</th><th>ClinVar</th><th>REVEL</th>" +
          scoreColumnsHead() +
          "</tr></thead><tbody>" +
          vs.map((v) => variantRow(v.gene || contextSymbol, v)).join("") +
          "</tbody></table></div>" +
          "</div>"
        );
      })
      .join("");

    return summary + sections;
  }

  function scoreColumnsHead() {
    return "<th>CADD</th><th>SpliceAI</th><th>AlphaMis.</th><th>QUAL</th><th>Reads (alt/total)</th>";
  }

  function scoreColumnsCells(v) {
    const caddFmt = v.cadd !== null && v.cadd !== undefined ? Number(v.cadd).toFixed(1) : "—";
    const caddHot = v.cadd && v.cadd >= 20;
    const caddBar = v.cadd !== null && v.cadd !== undefined ? Math.min(100, Math.round((v.cadd / 40) * 100)) : 0;

    const spliceFmt = v.spliceai !== null && v.spliceai !== undefined ? Number(v.spliceai).toFixed(2) : "—";
    const spliceHot = v.spliceai && v.spliceai >= 0.5;
    const spliceBar = v.spliceai !== null && v.spliceai !== undefined ? Math.min(100, Math.round(v.spliceai * 100)) : 0;

    const amFmt = v.alphamissense !== null && v.alphamissense !== undefined ? Number(v.alphamissense).toFixed(3) : "—";
    const amHot = v.alphamissense && v.alphamissense >= 0.564;
    const amBar = v.alphamissense !== null && v.alphamissense !== undefined ? Math.min(100, Math.round(v.alphamissense * 100)) : 0;

    const qualFmt = v.qual !== null && v.qual !== undefined ? Number(v.qual).toFixed(1) : "—";
    const readsFmt = v.reads && v.reads.total ? v.reads.matching + "/" + v.reads.total : "—";
    const readsPct = v.reads && v.reads.total ? Math.round((v.reads.matching / v.reads.total) * 100) + "%" : "";

    return (
      '<td><div class="score-cell"><div class="score-bar-track"><div class="score-bar-fill' + (caddHot ? " hot" : "") + '" style="width:' + caddBar + '%"></div></div><span class="mono">' + caddFmt + "</span></div></td>" +
      '<td><div class="score-cell"><div class="score-bar-track"><div class="score-bar-fill' + (spliceHot ? " hot" : "") + '" style="width:' + spliceBar + '%"></div></div><span class="mono">' + spliceFmt + "</span></div></td>" +
      '<td><div class="score-cell"><div class="score-bar-track"><div class="score-bar-fill' + (amHot ? " hot" : "") + '" style="width:' + amBar + '%"></div></div><span class="mono">' + amFmt + "</span></div></td>" +
      '<td class="mono">' + qualFmt + "</td>" +
      '<td><div class="score-cell"><span class="reads-frac">' + readsFmt + '</span><span class="reads-pct">' + readsPct + "</span></div></td>"
    );
  }

  function variantRow(geneSymbol, v) {
    const badgeCls = { concern: "badge--concern", protective: "badge--protect", uncertain: "badge--uncertain", uncategorized: "badge--neutral" };
    const rowKey = (geneSymbol || v.gene) + ":" + v.id;
    const isExpanded = state.expandedVariantRows.has(rowKey);

    const mainRow =
      '<tr data-gene="' + (geneSymbol || v.gene) + '" data-vid="' + v.id + '" class="' + (isExpanded ? "expanded" : "") + '">' +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);"><span class="caret" style="margin-right:6px;">' + (isExpanded ? "\u25BE" : "\u25B8") + '</span>' + v.id + "</td>" +
      '<td class="mono">' + (geneSymbol || v.gene) + "</td>" +
      '<td class="mono">' + v.genotype + "</td>" +
      "<td>" + v.zygosity + "</td>" +
      "<td>" + phaseTag(v.phase) + "</td>" +
      "<td>" + (typeof v.maf === 'number' ? v.maf.toFixed(4) : v.maf) + "</td>" +
      "<td>" + v.consequence.map((c) => '<span class="chip">' + c + "</span>").join("") + "</td>" +
      "<td><span class=\"badge " + (badgeCls[v.category] || "badge--neutral") + "\">" + v.clinvar + "</span></td>" +
      "<td>" + (v.revel === null || v.revel === undefined ? '<span class="na-cell">\u2014</span>' : v.revel) + "</td>" +
      scoreColumnsCells(v) +
      "</tr>";

    if (!isExpanded) return mainRow;

    // Expandable Study Row
    const detailRow =
      '<tr class="variant-detail-row"><td colspan="14">' +
      '<div class="study-card" style="margin:6px 0;background:#fff;">' +
        '<div class="study-card__head">' +
          '<div style="font-weight:700;color:var(--teal-dark);">Locus: ' + v.coordinate + ' · Call Support: ' + (v.reads ? v.reads.matching + "/" + v.reads.total + " reads" : "N/A") + '</div>' +
          '<div style="font-size:11.5px;color:var(--slate);">Last Evaluated: ' + (v.lastEvaluated || "2026-07-06") + '</div>' +
        '</div>' +
        (v.studies && v.studies.length
          ? v.studies.map((s) =>
              '<div style="padding:6px 0;border-top:1px dashed var(--line);margin-top:6px;">' +
                '<div style="font-size:12.5px;font-weight:600;">' + s.finding + '</div>' +
                '<div style="font-size:11.5px;color:var(--slate);">' + s.source + ' · ' + s.genotypeRelevance + '</div>' +
              '</div>'
            ).join("")
          : '<div style="font-size:12px;color:var(--slate);padding:4px 0;">No individual GWAS study correlations published for this exact SNP.</div>') +
      '</div>' +
      '</td></tr>';

    return mainRow + detailRow;
  }

  function renderStudiesTabForVariants(variants) {
    const studiesList = [];
    variants.forEach((v) => {
      (v.studies || []).forEach((s) => {
        studiesList.push({ variant: v, study: s });
      });
    });

    if (!studiesList.length) {
      return '<div class="pub-empty">No individual GWAS study records linked to the variants in this selection.</div>';
    }

    return '<div class="studies-list">' +
      studiesList.map(({ variant: v, study: s }) =>
        '<div class="study-card">' +
          '<div class="study-card__head">' +
            '<span class="study-card__variant" data-gene="' + (v.gene || '') + '" data-vid="' + v.id + '">' + v.id + ' (' + (v.gene || '') + ')</span>' +
            '<span class="study-card__condition">' + s.condition + '</span>' +
          '</div>' +
          '<div class="study-card__finding">' + s.finding + '</div>' +
          '<div class="study-card__foot">' +
            '<span>' + s.source + '</span>' +
            '<span>' + s.genotypeRelevance + '</span>' +
          '</div>' +
        '</div>'
      ).join("") +
      '</div>';
  }

  function wireVariantsPane(wrap) {
    wrap.querySelectorAll(".variant-table tbody tr[data-vid]").forEach((tr) => {
      tr.addEventListener("click", () => {
        const rowKey = tr.dataset.gene + ":" + tr.dataset.vid;
        if (state.expandedVariantRows.has(rowKey)) state.expandedVariantRows.delete(rowKey);
        else state.expandedVariantRows.add(rowKey);
        renderInspectionPanel();
      });
    });
  }

  // -----------------------------------------------------------------
  // Genome Browser Interactive Modal (SVG Exon & Lollipop Visualizer)
  // -----------------------------------------------------------------
  function openGenomeModal(symbol, variants) {
    const modal = byId("genome-modal");
    const mount = byId("genome-browser-mount");
    const locusEl = byId("genome-modal-locus");
    const detailEl = byId("genome-browser-variant-detail");

    const gene = geneBySymbol(symbol) || {};
    byId("genome-modal-title").textContent = "Genome Locus Track — " + symbol;
    locusEl.textContent = gene.chromosome || (variants[0] ? variants[0].coordinate : "chrX");
    detailEl.innerHTML = '<div class="variant-info-box">Click on any lollipop marker below to inspect variant amino acid consequence and call quality.</div>';

    // Build SVG Exon Track with Lollipops
    const width = 720, height = 180, padX = 40;
    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    svg.setAttribute("class", "genome-browser-svg");
    svg.setAttribute("width", width);
    svg.setAttribute("height", height);
    svg.setAttribute("viewBox", "0 0 " + width + " " + height);

    // Coordinate axis
    const axisY = 120;
    const axisLine = document.createElementNS(svgns, "line");
    axisLine.setAttribute("x1", padX); axisLine.setAttribute("y1", axisY);
    axisLine.setAttribute("x2", width - padX); axisLine.setAttribute("y2", axisY);
    axisLine.setAttribute("stroke", "var(--line-strong)");
    axisLine.setAttribute("stroke-width", "2");
    svg.appendChild(axisLine);

    // Mock Exons along the axis
    const exons = [
      { start: 0.05, end: 0.18 },
      { start: 0.26, end: 0.42 },
      { start: 0.50, end: 0.68 },
      { start: 0.76, end: 0.92 }
    ];
    exons.forEach((ex, i) => {
      const exX = padX + ex.start * (width - 2 * padX);
      const exW = (ex.end - ex.start) * (width - 2 * padX);
      const rect = document.createElementNS(svgns, "rect");
      rect.setAttribute("x", exX); rect.setAttribute("y", axisY - 10);
      rect.setAttribute("width", exW); rect.setAttribute("height", 20);
      rect.setAttribute("fill", "var(--teal-dim)");
      rect.setAttribute("stroke", "var(--teal)");
      rect.setAttribute("rx", "3");
      svg.appendChild(rect);

      const label = document.createElementNS(svgns, "text");
      label.setAttribute("x", exX + exW / 2); label.setAttribute("y", axisY + 4);
      label.setAttribute("text-anchor", "middle");
      label.setAttribute("font-size", "10px"); label.setAttribute("fill", "var(--teal-dark)");
      label.textContent = "Exon " + (i + 1);
      svg.appendChild(label);
    });

    // Variant Lollipops
    const varList = Array.isArray(variants) ? variants : [variants];
    const lollipopGroup = document.createElementNS(svgns, "g");
    lollipopGroup.setAttribute("class", "lollipops-group");

    varList.forEach((v, i) => {
      const relPos = varList.length === 1 ? 0.48 : (0.1 + (i / Math.max(1, varList.length - 1)) * 0.8);
      const posX = padX + relPos * (width - 2 * padX);
      const stickY = 40 + (i % 2) * 20;

      const stick = document.createElementNS(svgns, "line");
      stick.setAttribute("x1", posX); stick.setAttribute("y1", stickY);
      stick.setAttribute("x2", posX); stick.setAttribute("y2", axisY - 10);
      stick.setAttribute("stroke", "var(--slate-light)");
      stick.setAttribute("stroke-width", "1.5");
      stick.setAttribute("stroke-dasharray", "2,2");
      lollipopGroup.appendChild(stick);

      const colorMap = { concern: "var(--concern)", protective: "var(--protect)", uncertain: "var(--uncertain)", uncategorized: "var(--teal)" };
      const circle = document.createElementNS(svgns, "circle");
      circle.setAttribute("cx", posX); circle.setAttribute("cy", stickY);
      circle.setAttribute("r", 6);
      circle.setAttribute("fill", colorMap[v.category] || "var(--teal)");
      circle.setAttribute("stroke", "#fff");
      circle.setAttribute("stroke-width", "1.5");

      circle.addEventListener("click", () => {
        detailEl.innerHTML =
          '<div class="variant-info-box">' +
            '<div style="font-weight:700;color:var(--teal-dark);font-size:13px;">' + v.id + ' (' + v.genotype + ') · ' + v.coordinate + '</div>' +
            '<div style="font-size:12px;margin-top:3px;"><strong>ClinVar:</strong> ' + v.clinvar + ' | <strong>CADD:</strong> ' + (v.cadd || "N/A") + ' | <strong>Zygosity:</strong> ' + v.zygosity + ' (' + v.phase + ')</div>' +
            '<div style="font-size:11.5px;color:var(--slate);margin-top:2px;">Consequences: ' + v.consequence.join(", ") + '</div>' +
          '</div>';
      });

      lollipopGroup.appendChild(circle);
    });

    svg.appendChild(lollipopGroup);
    mount.innerHTML = "";
    mount.appendChild(svg);

    modal.classList.add("open");
  }

  // -----------------------------------------------------------------
  // Global Tables: Genes, Variants, Analysis, Reports
  // -----------------------------------------------------------------
  let genesSortKey = "symbol";
  let genesSortDir = 1;

  function renderGenesTable() {
    const searchEl = byId("genes-search");
    const q = searchEl ? searchEl.value.trim().toLowerCase() : "";
    let rows = GENES.map((g) => {
      const heteroCount = g.variants.filter((v) => v.zygosity === "Heterozygous").length;
      const phasedHetCount = g.variants.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length;
      const phasedPct = heteroCount ? Math.round((phasedHetCount / heteroCount) * 100) : 0;
      return {
        symbol: g.symbol,
        organSystem: g.organSystem,
        variantsDetected: g.variantsDetected,
        pathogenic: g.variants.filter((v) => v.category === "concern").length,
        phasedPct: phasedPct,
        hpoTermCount: g.hpoTermCount,
        goTermCount: g.goTermCount
      };
    }).filter((r) => !q || r.symbol.toLowerCase().includes(q) || r.organSystem.toLowerCase().includes(q));

    rows.sort((a, b) => {
      const av = a[genesSortKey], bv = b[genesSortKey];
      if (typeof av === "string") return av.localeCompare(bv) * genesSortDir;
      return (av - bv) * genesSortDir;
    });

    const tbody = byId("genes-table").querySelector("tbody");
    tbody.innerHTML = rows.map((r) =>
      '<tr data-symbol="' + r.symbol + '" style="cursor:pointer;">' +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);">' + r.symbol + "</td>" +
      "<td>" + r.organSystem + "</td>" +
      "<td>" + r.variantsDetected + "</td>" +
      "<td>" + (r.pathogenic ? '<span class="badge badge--concern">' + r.pathogenic + "</span>" : "0") + "</td>" +
      '<td><div class="bar-cell"><div class="bar-track"><div class="bar-fill" style="width:' + r.phasedPct + '%"></div></div>' + r.phasedPct + "%</div></td>" +
      "<td>" + r.hpoTermCount + "</td>" +
      "<td>" + r.goTermCount + "</td>" +
      "</tr>"
    ).join("");

    tbody.querySelectorAll("tr").forEach((tr) => {
      tr.addEventListener("click", () => {
        const sym = tr.dataset.symbol;
        document.querySelectorAll(".tabs button").forEach((b) => b.classList.remove("active"));
        document.querySelector('.tabs button[data-view="ontology"]').classList.add("active");
        document.querySelectorAll(".view").forEach((v) => v.classList.remove("active"));
        byId("view-ontology").classList.add("active");
        selectGene(sym);
      });
    });
  }

  function allVariantsFlat() {
    return GENES.flatMap((g) => g.variants.map((v) => Object.assign({}, v, { gene: g.symbol })));
  }

  function renderVariantsTable() {
    const geneF = byId("filter-gene").value;
    const clinvarF = byId("filter-clinvar").value;
    const phaseF = byId("filter-phase").value;
    const zygF = byId("filter-zygosity").value;

    // Populate gene filter if empty
    const geneSelect = byId("filter-gene");
    if (geneSelect && geneSelect.children.length <= 1) {
      GENES.forEach((g) => {
        const opt = document.createElement("option");
        opt.value = g.symbol;
        opt.textContent = g.symbol;
        geneSelect.appendChild(opt);
      });
    }

    const rows = allVariantsFlat().filter((v) =>
      (!geneF || v.gene === geneF) &&
      (!clinvarF || v.clinvar === clinvarF) &&
      (!phaseF || v.phase === phaseF) &&
      (!zygF || v.zygosity === zygF)
    );

    const badgeCls = { concern: "badge--concern", protective: "badge--protect", uncertain: "badge--uncertain", uncategorized: "badge--neutral" };

    const tbody = byId("variants-table").querySelector("tbody");
    tbody.innerHTML = rows.map((v) =>
      '<tr data-gene="' + v.gene + '" data-vid="' + v.id + '" style="cursor:pointer;">' +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);">' + v.id + "</td>" +
      '<td class="mono">' + v.gene + "</td>" +
      "<td><span class=\"badge " + (badgeCls[v.category] || "badge--neutral") + "\">" + v.clinvar + "</span></td>" +
      "<td>" + (v.revel === null || v.revel === undefined ? '<span class="na-cell">\u2014</span>' : v.revel) + "</td>" +
      scoreColumnsCells(v) +
      '<td class="mono">' + v.coordinate + "</td>" +
      "<td>" + v.zygosity + "</td>" +
      "<td>" + phaseTag(v.phase) + "</td>" +
      "<td>" + (typeof v.maf === 'number' ? v.maf.toFixed(4) : v.maf) + "</td>" +
      "<td>" + v.consequence.map((c) => '<span class="chip">' + c + "</span>").join("") + "</td>" +
      "</tr>"
    ).join("");

    tbody.querySelectorAll("tr").forEach((tr) => {
      tr.addEventListener("click", () => {
        const sym = tr.dataset.gene;
        document.querySelectorAll(".tabs button").forEach((b) => b.classList.remove("active"));
        document.querySelector('.tabs button[data-view="ontology"]').classList.add("active");
        document.querySelectorAll(".view").forEach((v) => v.classList.remove("active"));
        byId("view-ontology").classList.add("active");
        selectGene(sym);
      });
    });
  }

  function renderAnalysisView() {
    const allVars = allVariantsFlat();
    const totalPath = allVars.filter((v) => v.category === "concern").length;
    const totalHetero = allVars.filter((v) => v.zygosity === "Heterozygous").length;
    const totalPhased = allVars.filter((v) => v.zygosity === "Heterozygous" && v.phase && v.phase !== "Unknown").length;
    const phasedPct = totalHetero ? Math.round((totalPhased / totalHetero) * 100) : 0;

    byId("analysis-kpis").innerHTML =
      kpi(allVars.length.toLocaleString(), "Total Actionable Variants") +
      kpi(totalPath, "Pathogenic / LP calls") +
      kpi(phasedPct + "%", "Heterozygous calls phased") +
      kpi(GENES.length, "Genes in panel");

    byId("prs-grid").innerHTML = PRS.map((p) =>
      '<div class="prs-card"><div class="prs-card__head"><span class="trait">' + p.trait + '</span><span class="badge badge--' +
      (p.category === "HIGH" ? "concern" : p.category === "PROTECTIVE" ? "protect" : "uncertain") + '">' + p.category + '</span></div>' +
      '<div class="prs-track"><div class="prs-fill ' + p.category + '" style="width:' + p.percentile + '%"></div></div>' +
      '<div class="prs-foot"><span>' + p.percentile + "th percentile · " + p.organSystem + '</span><a href="#">' + p.pgsId + '</a></div>' +
      "</div>"
    ).join("");

    byId("pgx-body").innerHTML = PGX.map((p) =>
      "<tr><td class=\"mono\" style=\"font-weight:700;\">" + p.gene + '</td><td class="mono">' + p.diplotype + "</td><td>" + p.phenotype +
      "</td><td>" + p.drug + '</td><td><span class="action-tier ' + p.actionTier + '">' + p.actionTier + "</span></td><td>" + p.recommendation + "</td></tr>"
    ).join("");
  }

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
  // Helper Formatters
  // -----------------------------------------------------------------
  function kpi(value, label) {
    return '<div class="kpi"><div class="val">' + esc(value) + '</div><div class="lbl">' + esc(label) + '</div></div>';
  }

  function phaseTag(phase) {
    if (!phase || phase === "Unknown") return '<span class="badge badge--phase-unknown">Unknown</span>';
    if (phase === "Maternal") return '<span class="badge badge--phase-mat">Maternal</span>';
    if (phase === "Paternal") return '<span class="badge badge--phase-pat">Paternal</span>';
    return '<span class="badge badge--neutral">' + phase + '</span>';
  }

  function pubCard(p) {
    return (
      '<div class="pub-card">' +
      '<div class="pub-card__title">' + p.title + '</div>' +
      '<div class="pub-card__authors">' + p.authors + ' · <span class="mono">' + p.journal + ' (' + p.year + ')</span></div>' +
      '<div class="pub-card__relevance">' + p.relevance + '</div>' +
      '<a class="pub-card__link" href="https://pubmed.ncbi.nlm.nih.gov/' + p.pmid + '/" target="_blank" rel="noopener">PMID ' + p.pmid + ' ↗</a>' +
      '</div>'
    );
  }

  // Close genome browser modal
  const modalClose = byId("genome-modal-close");
  if (modalClose) {
    modalClose.addEventListener("click", () => {
      byId("genome-modal").classList.remove("open");
    });
  }

  // -----------------------------------------------------------------
  // Bootstrapping
  // -----------------------------------------------------------------
  document.addEventListener("DOMContentLoaded", () => {
    initTabs();
    initOntologySwitch();

    // Setup filter listeners in variants view
    ["filter-gene", "filter-clinvar", "filter-phase", "filter-zygosity"].forEach((id) => {
      const el = byId(id);
      if (el) el.addEventListener("change", renderVariantsTable);
    });

    const geneSearch = byId("genes-search");
    if (geneSearch) geneSearch.addEventListener("input", renderGenesTable);

    // Initial render
    renderLeftPanel();

    // Default select first Level 1 Category
    const filtered = filteredOntology();
    if (filtered.length > 0) {
      selectCategory({
        type: "category",
        level: 1,
        id: filtered[0].group.id,
        label: filtered[0].group.label,
        levelName: "Organ System / Root Category",
        genes: filtered[0].genes,
        subcategories: filtered[0].subcategories
      });
    }
  });

})();
