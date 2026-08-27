/**
 * Genomic Ontology Explorer — app logic (v4.5)
 * Gene-grouped Tiered Variants with NCBI descriptions,
 * Non-redundant Comprehensive Variant Details (SpliceAI AG/AL/DG/DL, AlphaMissense, ACMG, cDNA/Protein),
 * Interactive Scalable Graph with Zoom Controls, and Persistent Subtab Navigation.
 */

(function () {
  "use strict";

  // -----------------------------------------------------------------
  // Shared state
  // -----------------------------------------------------------------
  const state = {
    tree: "hpo",
    layout: "tree",               // default to Tree for rich multi-level exploration
    scopeFindingsOnly: false,
    treeSearch: "",
    selectedTarget: null,         // node or gene object
    activeSubTab: "overview",     // maintains active subtab (overview, phenotypes, variants, studies, publications)
    graphScale: 1.0,              // graph zoom scale
    expandedNodes: new Set(),      // ontology node ids currently expanded
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
        
        const groups = filteredOntology();
        if (groups && groups.length) {
          selectCategory(groups[0]);
        }
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
  // Filtered Ontology Data Model (Recursive Arbitrary-Depth)
  // -----------------------------------------------------------------
  function filterNode(node) {
    let nodeGenes = Array.from(new Set(node.genes || []));
    if (state.scopeFindingsOnly) nodeGenes = nodeGenes.filter(geneHasFindings);

    const filteredChildren = (node.children || [])
      .map(filterNode)
      .filter(Boolean);

    const nodeHit = matchesSearch(node.label) || matchesSearch(node.id) || nodeGenes.some(matchesSearch);
    if (state.scopeFindingsOnly && !nodeGenes.length && !filteredChildren.length) return null;
    if (state.treeSearch && !nodeHit && !filteredChildren.length) return null;

    return {
      ...node,
      genes: nodeGenes,
      children: filteredChildren
    };
  }

  function filteredOntology() {
    const ont = ONTOLOGIES[state.tree] || ONTOLOGIES["hpo"];
    if (!ont || !ont.groups) return [];
    return ont.groups.map(filterNode).filter(Boolean);
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
  // Recursive Multi-Level Tree & List Rendering
  // -----------------------------------------------------------------
  function renderTreeOrList() {
    const data = filteredOntology();
    const scroll = byId("tree-scroll");
    scroll.className = "tree-scroll" + (state.layout === "list" ? " tree-list-mode" : "");
    scroll.innerHTML = "";

    if (state.layout === "list") {
      data.forEach((group) => {
        const groupEl = document.createElement("div");
        groupEl.className = "tree-group";
        const isOpen = state.expandedNodes.has(group.id) || !!state.treeSearch;
        const isSelected = state.selectedTarget && state.selectedTarget.id === group.id;

        const row = document.createElement("div");
        row.className = "tree-row tree-row--organ" + (isOpen ? " expanded" : "") + (isSelected ? " selected" : "");
        row.innerHTML =
          '<span class="caret">' + (isOpen ? "\u25BE" : "\u25B8") + '</span><span class="node-dot"></span>' +
          "<span>" + group.label + "</span>" +
          '<span class="count-chip">' + group.genes.length + " genes</span>";

        row.addEventListener("click", () => {
          if (state.expandedNodes.has(group.id)) state.expandedNodes.delete(group.id);
          else state.expandedNodes.add(group.id);
          selectCategory(group);
        });
        groupEl.appendChild(row);

        const childWrap = document.createElement("div");
        childWrap.className = "tree-children" + (isOpen ? " open" : "");
        group.genes.forEach((sym) => childWrap.appendChild(makeGeneRow(sym)));
        groupEl.appendChild(childWrap);
        scroll.appendChild(groupEl);
      });
    } else {
      data.forEach((rootNode) => {
        scroll.appendChild(renderTreeNodeRecursive(rootNode, 0));
      });
    }

    if (!scroll.children.length) {
      scroll.innerHTML = '<div style="padding:20px;color:var(--slate);font-size:12.5px;">No matches under the current filters.</div>';
    }
  }

  function renderTreeNodeRecursive(node, depth) {
    const wrap = document.createElement("div");
    wrap.className = "tree-node-wrap";

    const isOpen = state.expandedNodes.has(node.id) || (depth === 0 && state.expandedNodes.size === 0) || !!state.treeSearch;
    const isSelected = state.selectedTarget && state.selectedTarget.id === node.id;
    const hasChildren = (node.children && node.children.length > 0) || (node.genes && node.genes.length > 0);

    const row = document.createElement("div");
    const levelClass = depth === 0 ? "tree-row--organ" : (depth === 1 ? "tree-row--subcat" : "tree-row--term");
    row.className = "tree-row " + levelClass + (isOpen ? " expanded" : "") + (isSelected ? " selected" : "");
    row.style.paddingLeft = (12 + depth * 14) + "px";

    const dotStyle = depth === 0 ? "" : (depth === 1 ? "background:var(--teal);" : "background:var(--slate-light);");
    const countLabel = node.genes.length + (depth === 0 ? " genes" : "");

    row.innerHTML =
      (hasChildren ? '<span class="caret">' + (isOpen ? "\u25BE" : "\u25B8") + '</span>' : '<span class="caret" style="opacity:0;">•</span>') +
      '<span class="node-dot" style="' + dotStyle + '"></span>' +
      '<span class="tree-node-title">' + node.label + '</span>' +
      (node.id && (node.id.startsWith("HP:") || node.id.startsWith("GO:") || node.id.startsWith("ORGAN:")) ? '<span class="count-chip mono">' + node.id + '</span>' : '') +
      '<span class="count-chip">' + countLabel + '</span>';

    row.addEventListener("click", (e) => {
      e.stopPropagation();
      if (hasChildren) {
        if (state.expandedNodes.has(node.id)) state.expandedNodes.delete(node.id);
        else state.expandedNodes.add(node.id);
      }
      selectCategory(node);
    });

    wrap.appendChild(row);

    const childWrap = document.createElement("div");
    childWrap.className = "tree-children" + (isOpen ? " open" : "");

    if (node.children && node.children.length > 0) {
      node.children.forEach((child) => {
        childWrap.appendChild(renderTreeNodeRecursive(child, depth + 1));
      });
    } else if (node.genes && node.genes.length > 0) {
      node.genes.forEach((sym) => {
        const geneRow = makeGeneRow(sym);
        geneRow.style.paddingLeft = (16 + (depth + 1) * 14) + "px";
        childWrap.appendChild(geneRow);
      });
    }

    wrap.appendChild(childWrap);
    return wrap;
  }

  function makeGeneRow(symbol) {
    const g = geneBySymbol(symbol);
    const row = document.createElement("div");
    const isSelected = state.selectedTarget && state.selectedTarget.type === "gene" && state.selectedTarget.symbol === symbol;
    row.className = "tree-row tree-row--gene" + (isSelected ? " selected" : "");
    const flagged = g && g.variants.some((v) => v.category === "concern");
    row.innerHTML =
      '<span class="node-dot"></span><span style="font-weight:600;">' + symbol + "</span>" +
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

  function selectCategory(node) {
    state.selectedTarget = { ...node, type: "category" };
    if (state.layout === "graph") renderGraph(); else renderTreeOrList();
    renderInspectionPanel();
  }

  // -----------------------------------------------------------------
  // Graph layout mode — Scalable, Zoomable SVG Hierarchy
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

    // Toolbar for Zoom In / Zoom Out / Reset
    const toolbar = document.createElement("div");
    toolbar.className = "graph-toolbar";
    toolbar.innerHTML =
      '<button id="graph-zoom-in" title="Zoom In">+</button>' +
      '<button id="graph-zoom-out" title="Zoom Out">−</button>' +
      '<button id="graph-zoom-reset" title="Reset View">↺</button>';
    host.appendChild(toolbar);

    toolbar.querySelector("#graph-zoom-in").addEventListener("click", (e) => {
      e.stopPropagation();
      state.graphScale = Math.min(2.0, state.graphScale + 0.15);
      applyGraphTransform();
    });
    toolbar.querySelector("#graph-zoom-out").addEventListener("click", (e) => {
      e.stopPropagation();
      state.graphScale = Math.max(0.4, state.graphScale - 0.15);
      applyGraphTransform();
    });
    toolbar.querySelector("#graph-zoom-reset").addEventListener("click", (e) => {
      e.stopPropagation();
      state.graphScale = 1.0;
      applyGraphTransform();
    });

    const data = filteredOntology();
    const ROW_H = 34, L1_X = 30, L2_X = 190, L3_X = 355, GENE_X = 520, PAD_TOP = 25, BASE_WIDTH = 800;

    const canvasWrap = document.createElement("div");
    canvasWrap.id = "graph-canvas-wrap";
    canvasWrap.style.transformOrigin = "top left";
    canvasWrap.style.transition = "transform 0.1s ease-out";

    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    svg.setAttribute("class", "graph-svg");
    svg.setAttribute("width", BASE_WIDTH);

    const geneChips = [];
    let y = PAD_TOP;
    const edges = [];
    const nodes = [];

    data.forEach((root) => {
      const rootY = y;
      nodes.push({ id: root.id, x: L1_X, y: rootY, label: root.label, node: root, level: 1 });
      y += ROW_H;

      (root.children || []).forEach((sub) => {
        const subY = y;
        nodes.push({ id: sub.id, x: L2_X, y: subY, label: sub.label, node: sub, level: 2 });
        edges.push({ x1: L1_X + 8, y1: rootY, x2: L2_X - 8, y2: subY, strokeWidth: "1.75" });
        y += ROW_H;

        (sub.children || []).slice(0, 3).forEach((term) => {
          const termY = y;
          nodes.push({ id: term.id, x: L3_X, y: termY, label: term.label, node: term, level: 3 });
          edges.push({ x1: L2_X + 8, y1: subY, x2: L3_X - 8, y2: termY, strokeWidth: "1.5" });
          geneChips.push({ y: termY, genes: term.genes || [] });
          y += ROW_H;
        });
      });
      y += 14;
    });

    const totalHeight = Math.max(500, y + 40);
    svg.setAttribute("height", totalHeight);
    svg.setAttribute("viewBox", "0 0 " + BASE_WIDTH + " " + totalHeight);

    // Balanced Contrast Lines (#64748b with subtle curvature)
    edges.forEach((e) => {
      const path = document.createElementNS(svgns, "path");
      const midX = (e.x1 + e.x2) / 2;
      const d = `M ${e.x1} ${e.y1} C ${midX} ${e.y1}, ${midX} ${e.y2}, ${e.x2} ${e.y2}`;
      path.setAttribute("d", d);
      path.setAttribute("fill", "none");
      path.setAttribute("stroke", "#64748b"); // balanced slate-500 line
      path.setAttribute("stroke-width", e.strokeWidth || "1.75");
      path.setAttribute("stroke-opacity", "0.75");
      svg.appendChild(path);
    });

    // Draw Nodes
    nodes.forEach((n) => {
      const isSel = state.selectedTarget && state.selectedTarget.id === n.id;
      const gNode = document.createElementNS(svgns, "g");
      gNode.style.cursor = "pointer";

      const circle = document.createElementNS(svgns, "circle");
      circle.setAttribute("cx", n.x); circle.setAttribute("cy", n.y);
      circle.setAttribute("r", n.level === 1 ? 6.5 : (n.level === 2 ? 5.0 : 4.0));
      circle.setAttribute("fill", isSel ? "var(--concern)" : (n.level === 1 ? "var(--teal-dark)" : (n.level === 2 ? "var(--teal)" : "#64748b")));
      circle.setAttribute("stroke", "#ffffff");
      circle.setAttribute("stroke-width", "1.5");
      gNode.appendChild(circle);

      const text = document.createElementNS(svgns, "text");
      text.setAttribute("x", n.x + 12);
      text.setAttribute("y", n.y + 4);
      text.setAttribute("font-size", n.level === 1 ? "12px" : (n.level === 2 ? "11px" : "10px"));
      text.setAttribute("font-weight", n.level === 1 ? "700" : "600");
      text.setAttribute("fill", isSel ? "var(--teal-dark)" : "var(--ink)");
      text.textContent = n.label.length > 24 ? n.label.substring(0, 22) + "…" : n.label;
      gNode.appendChild(text);

      gNode.addEventListener("click", () => selectCategory(n.node));
      svg.appendChild(gNode);
    });

    canvasWrap.appendChild(svg);

    // Gene chips layer
    const geneLayer = document.createElement("div");
    geneLayer.className = "graph-genes-layer";
    geneLayer.style.width = BASE_WIDTH + "px";
    geneLayer.style.height = totalHeight + "px";
    geneChips.forEach((row) => {
      (row.genes || []).slice(0, 4).forEach((sym, i) => {
        const chip = document.createElement("div");
        const isGeneSel = state.selectedTarget && state.selectedTarget.type === "gene" && state.selectedTarget.symbol === sym;
        chip.className = "graph-gene-chip" + (isGeneSel ? " selected" : "");
        chip.textContent = sym;
        chip.style.left = (GENE_X + i * 66) + "px";
        chip.style.top = (row.y - 10) + "px";
        chip.addEventListener("click", () => selectGene(sym));
        geneLayer.appendChild(chip);
      });
    });
    canvasWrap.appendChild(geneLayer);

    host.appendChild(canvasWrap);
    applyGraphTransform();
  }

  function applyGraphTransform() {
    const wrap = byId("graph-canvas-wrap");
    if (wrap) {
      wrap.style.transform = "scale(" + state.graphScale + ")";
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

    const activeTab = state.activeSubTab || "overview";

    wrap.innerHTML =
      '<div class="gene-head">' +
        '<div>' +
          '<h1><span class="sym">' + g.symbol + '</span> <span class="full">' + g.name + '</span></h1>' +
          '<div class="gene-coord mono" style="font-size:11.5px;color:var(--slate);margin-top:2px;">' + g.chromosome + " · pLI " + g.pli + " · LOEUF " + g.loeuf + " · " + g.organSystem + '</div>' +
        '</div>' +
        '<button class="btn btn--primary" id="genome-browser-btn">View in Genome Browser</button>' +
      '</div>' +
      '<div class="gene-tabs" id="gene-tabs">' +
        '<button data-tab="overview" class="' + (activeTab === "overview" ? "active" : "") + '">Overview</button>' +
        '<button data-tab="phenotypes" class="' + (activeTab === "phenotypes" ? "active" : "") + '">Phenotypes (' + g.hpoTermCount + ')</button>' +
        '<button data-tab="variants" class="' + (activeTab === "variants" ? "active" : "") + '">Variants (' + g.variants.length + ')</button>' +
        '<button data-tab="studies" class="' + (activeTab === "studies" ? "active" : "") + '">Studies</button>' +
        '<button data-tab="publications" class="' + (activeTab === "publications" ? "active" : "") + '">Publications (' + (g.publications ? g.publications.length : 0) + ')</button>' +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "overview" ? "active" : "") + '" data-pane="overview">' +
        '<div class="kpi-row">' +
          kpi(g.variantsDetected, "Variants detected") +
          kpi(pathCount, "Potential concerns") +
          kpi(protectCount, "Protective factors") +
          kpi(uncertainCount, "Uncertain") +
          kpi(phasedPct + "%", "Het. variants phased") +
          kpi(g.hpoTermCount, "HPO terms") +
          kpi(g.goTermCount, "GO terms") +
          kpi(g.publications ? g.publications.length : 0, "Publications") +
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

      '<div class="gene-pane ' + (activeTab === "phenotypes" ? "active" : "") + '" data-pane="phenotypes">' +
        (g.hpoTerms && g.hpoTerms.length
          ? '<div class="hpo-card-grid">' + g.hpoTerms.map((t) =>
              '<div class="hpo-card"><div class="id">' + t.id + '</div><div class="label">' + t.label + '</div><div class="evidence">' + t.evidence + '</div></div>'
            ).join("") + '</div>'
          : '<div class="pub-empty">No curated HPO associations recorded for this gene.</div>') +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "variants" ? "active" : "") + '" data-pane="variants">' + renderVariantSectionList(g.variants, g.symbol) + '</div>' +

      '<div class="gene-pane ' + (activeTab === "studies" ? "active" : "") + '" data-pane="studies">' + renderStudiesTabForVariants(g.variants) + '</div>' +

      '<div class="gene-pane ' + (activeTab === "publications" ? "active" : "") + '" data-pane="publications">' +
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
  function collectDescendantSubOntologies(node) {
    const subTerms = [];
    function walk(n) {
      if (n.children && n.children.length) {
        n.children.forEach((c) => {
          subTerms.push(c);
          walk(c);
        });
      }
    }
    walk(node);
    return subTerms;
  }

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

    const subOntologies = collectDescendantSubOntologies(cat);
    const allPubs = geneObjects.flatMap((g) => g.publications || []);

    const levelBadge = cat.level ? "Level " + cat.level : (cat.id && (cat.id.startsWith("HP:") || cat.id.startsWith("GO:") || cat.id.startsWith("ORGAN:")) ? cat.id : "Category");

    const activeTab = state.activeSubTab || "overview";

    wrap.innerHTML =
      '<div class="category-head">' +
        '<div>' +
          '<h1><span>' + cat.label + '</span> <span class="category-level-badge">' + levelBadge + '</span></h1>' +
          '<div class="category-meta mono">' + (cat.id ? cat.id + " · " : "") + geneObjects.length + " Genes · " + totalVariants + " Actionable Variants</div>" +
        '</div>' +
      '</div>' +

      '<div class="gene-tabs" id="gene-tabs">' +
        '<button data-tab="overview" class="' + (activeTab === "overview" ? "active" : "") + '">Overview</button>' +
        '<button data-tab="phenotypes" class="' + (activeTab === "phenotypes" ? "active" : "") + '">Phenotypes / Sub-Ontologies (' + subOntologies.length + ')</button>' +
        '<button data-tab="variants" class="' + (activeTab === "variants" ? "active" : "") + '">Variants (' + totalVariants + ')</button>' +
        '<button data-tab="studies" class="' + (activeTab === "studies" ? "active" : "") + '">Studies</button>' +
        '<button data-tab="publications" class="' + (activeTab === "publications" ? "active" : "") + '">Publications (' + allPubs.length + ')</button>' +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "overview" ? "active" : "") + '" data-pane="overview">' +
        '<div class="kpi-row">' +
          kpi(geneObjects.length, "Genes in domain") +
          kpi(pathCount, "Potential concerns") +
          kpi(protectCount, "Protective factors") +
          kpi(uncertainCount, "Uncertain") +
          kpi(phasedPct + "%", "Het. variants phased") +
          kpi(subOntologies.length, "Sub-ontologies") +
          kpi(totalVariants, "Total variants") +
          kpi(allPubs.length, "Publications") +
        '</div>' +

        '<div class="panel-box">' +
          '<h3>Clinical & Biological Domain Summary</h3>' +
          '<p>Roll-up analysis of findings categorized under <strong>' + cat.label + '</strong> (' + (cat.id || "Domain") + '). Encompasses <strong>' + geneObjects.length + ' genes</strong> with <strong>' + totalVariants + ' actionable variants</strong> identified across patient sequencing data.</p>' +
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

      '<div class="gene-pane ' + (activeTab === "phenotypes" ? "active" : "") + '" data-pane="phenotypes">' +
        '<div class="section-title"><span class="eyebrow">Sub-Level Ontologies & Phenotypic Branches (' + subOntologies.length + ')</span></div>' +
        (subOntologies.length
          ? '<div class="hpo-card-grid">' + subOntologies.map((s) =>
              '<div class="hpo-card" style="cursor:pointer;" data-node-id="' + s.id + '">' +
                '<div class="id" style="display:flex;justify-content:space-between;"><span>' + s.id + '</span><span class="category-level-badge" style="font-size:9px;">Level ' + (s.level || "Sub") + '</span></div>' +
                '<div class="label" style="font-weight:700;margin-top:4px;">' + s.label + '</div>' +
                '<div class="evidence" style="margin-top:6px;font-size:11px;color:var(--teal-dark);font-weight:600;">' + (s.genes ? s.genes.length + ' Genes: ' + s.genes.slice(0, 8).join(", ") + (s.genes.length > 8 ? "…" : "") : "") + '</div>' +
              '</div>'
            ).join("") + '</div>'
          : '<div class="pub-empty">This is a terminal leaf phenotype term. View associated genes above or in the Variants tab.</div>') +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "variants" ? "active" : "") + '" data-pane="variants">' +
        renderVariantSectionList(allVariants, cat.label) +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "studies" ? "active" : "") + '" data-pane="studies">' +
        renderStudiesTabForVariants(allVariants) +
      '</div>' +

      '<div class="gene-pane ' + (activeTab === "publications" ? "active" : "") + '" data-pane="publications">' +
        '<div class="pub-grid">' +
          allPubs.map(pubCard).join("") +
        '</div>' +
        (!allPubs.length ? '<div class="pub-empty">No curated publications directly indexed for this category.</div>' : "") +
      '</div>';

    wireTabSwitching(wrap);

    wrap.querySelectorAll(".gene-chip-item").forEach((chip) => {
      chip.addEventListener("click", () => selectGene(chip.dataset.gene));
    });

    wrap.querySelectorAll(".hpo-card[data-node-id]").forEach((card) => {
      card.addEventListener("click", () => {
        const targetNodeId = card.dataset.nodeId;
        const targetSub = subOntologies.find((s) => s.id === targetNodeId);
        if (targetSub) {
          state.expandedNodes.add(targetNodeId);
          selectCategory(targetSub);
        }
      });
    });

    wireVariantsPane(wrap);
  }

  function wireTabSwitching(wrap) {
    wrap.querySelectorAll("#gene-tabs button").forEach((btn) => {
      btn.addEventListener("click", () => {
        state.activeSubTab = btn.dataset.tab; // Store active subtab to prevent jump
        wrap.querySelectorAll("#gene-tabs button").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        wrap.querySelectorAll(".gene-pane").forEach((p) => p.classList.remove("active"));
        const targetPane = wrap.querySelector('.gene-pane[data-pane="' + btn.dataset.tab + '"]');
        if (targetPane) targetPane.classList.add("active");
      });
    });
  }

  // -----------------------------------------------------------------
  // Tiered Variants List: Grouped by Gene with NCBI/OMIM Descriptions
  // -----------------------------------------------------------------
  const CATEGORY_META = {
    concern: { label: "Potential concerns", cls: "concern" },
    protective: { label: "Protective associations", cls: "protect" },
    uncertain: { label: "Uncertain significance", cls: "uncertain" },
    uncategorized: { label: "Not categorized / Baseline", cls: "neutral" }
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

        // Group variants by Gene within this category
        const byGene = {};
        vs.forEach((v) => {
          const sym = v.gene || contextSymbol;
          if (!byGene[sym]) byGene[sym] = [];
          byGene[sym].push(v);
        });

        const geneBlocks = Object.keys(byGene).map((sym) => {
          const g = geneBySymbol(sym) || {};
          const geneVars = byGene[sym];
          return (
            '<div class="variant-gene-group">' +
              '<div class="variant-gene-header">' +
                '<div>' +
                  '<div>' +
                    '<span class="variant-gene-title">' + sym + '</span>' +
                    '<span class="variant-gene-fullname">' + (g.name || '') + '</span>' +
                    '<span class="count-chip" style="margin-left:8px;">' + geneVars.length + ' variant' + (geneVars.length === 1 ? '' : 's') + '</span>' +
                  '</div>' +
                  '<div class="variant-gene-desc">' + (g.summary || 'Essential clinical genomic locus.') + '</div>' +
                '</div>' +
                '<div style="flex-shrink:0;">' +
                  '<a class="ref-link" href="' + (g.links ? g.links.ncbiGene : 'https://www.ncbi.nlm.nih.gov/gene/?term=' + sym) + '" target="_blank" rel="noopener" style="font-size:11px;padding:3px 8px;">NCBI ↗</a>' +
                '</div>' +
              '</div>' +
              '<div class="data-table-wrap">' +
                '<table class="variant-table">' +
                  '<thead><tr>' +
                    '<th>Variant ID</th><th>Genotype</th><th>Zygosity</th><th>Phase</th><th>MAF</th><th>Consequence</th><th>ClinVar</th><th>REVEL</th>' +
                    scoreColumnsHead() +
                  '</tr></thead>' +
                  '<tbody>' +
                    geneVars.map((v) => variantRow(sym, v)).join("") +
                  '</tbody>' +
                '</table>' +
              '</div>' +
            '</div>'
          );
        }).join("");

        return (
          '<div class="variant-section" style="margin-bottom:24px;">' +
            '<div class="variant-section__head ' + meta.cls + '"><span>' + vs.length + " " + meta.label.toUpperCase() + '</span><span class="caret">\u25BE</span></div>' +
            geneBlocks +
          '</div>'
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

  // -----------------------------------------------------------------
  // Non-Redundant Comprehensive Variant Detail Drawer
  // -----------------------------------------------------------------
  function variantRow(geneSymbol, v) {
    const badgeCls = { concern: "badge--concern", protective: "badge--protect", uncertain: "badge--uncertain", uncategorized: "badge--neutral" };
    const rowKey = (geneSymbol || v.gene) + ":" + v.id;
    const isExpanded = state.expandedVariantRows.has(rowKey);

    const mainRow =
      '<tr data-gene="' + (geneSymbol || v.gene) + '" data-vid="' + v.id + '" class="' + (isExpanded ? "expanded" : "") + '">' +
      '<td class="mono" style="font-weight:700;color:var(--teal-dark);"><span class="caret" style="margin-right:6px;">' + (isExpanded ? "\u25BE" : "\u25B8") + '</span>' + v.id + "</td>" +
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

    // SpliceAI Delta Breakdown
    const sDetails = v.spliceaiDetails || {};
    const agStr = sDetails.ag !== undefined && sDetails.ag !== null ? Number(sDetails.ag).toFixed(2) : "0.00";
    const alStr = sDetails.al !== undefined && sDetails.al !== null ? Number(sDetails.al).toFixed(2) : "0.00";
    const dgStr = sDetails.dg !== undefined && sDetails.dg !== null ? Number(sDetails.dg).toFixed(2) : "0.00";
    const dlStr = sDetails.dl !== undefined && sDetails.dl !== null ? Number(sDetails.dl).toFixed(2) : "0.00";

    const detailRow =
      '<tr class="variant-detail-row"><td colspan="13">' +
      '<div class="variant-drawer-grid">' +
        // Box 1: Genomic & Transcript Context
        '<div class="variant-drawer-box">' +
          '<h4>Genomic & Transcript</h4>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Coordinates:</span><span class="drawer-metric-val">' + v.coordinate + '</span></div>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Transcript:</span><span class="drawer-metric-val">' + (v.transcript || 'Canonical') + '</span></div>' +
          (v.cchange ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">cDNA Change:</span><span class="drawer-metric-val">' + v.cchange + '</span></div>' : '') +
          (v.achange ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Protein Change:</span><span class="drawer-metric-val">' + v.achange + '</span></div>' : '') +
          (v.vaf !== undefined && v.vaf !== null ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Allele Fraction (VAF):</span><span class="drawer-metric-val">' + Number(v.vaf).toFixed(2) + '</span></div>' : '') +
        '</div>' +

        // Box 2: Functional In-Silico Scores
        '<div class="variant-drawer-box">' +
          '<h4>In-Silico Predictions</h4>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">AlphaMissense:</span><span class="drawer-metric-val">' + (v.alphamissense !== null && v.alphamissense !== undefined ? Number(v.alphamissense).toFixed(3) + ' (' + (v.amClass || 'evaluated') + ')' : '—') + '</span></div>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">CADD Phred:</span><span class="drawer-metric-val">' + (v.cadd ? Number(v.cadd).toFixed(1) : '—') + '</span></div>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">REVEL Score:</span><span class="drawer-metric-val">' + (v.revel !== null && v.revel !== undefined ? v.revel : '—') + '</span></div>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">SpliceAI Delta:</span><span class="drawer-metric-val">AG:' + agStr + ' AL:' + alStr + ' DG:' + dgStr + ' DL:' + dlStr + '</span></div>' +
        '</div>' +

        // Box 3: ACMG & ClinVar Evidence
        '<div class="variant-drawer-box">' +
          '<h4>ClinVar & ACMG Criteria</h4>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Classification:</span><span class="drawer-metric-val">' + v.clinvar + '</span></div>' +
          '<div class="drawer-metric-row"><span class="drawer-metric-lbl">Review Status:</span><span class="drawer-metric-val">' + (v.clinvarRev || 'criteria provided') + '</span></div>' +
          (v.acmgPm5 ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PM5:</span><span class="drawer-metric-val">Known pathogenic missense at codon</span></div>' : '') +
          (v.acmgPs1 ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PS1:</span><span class="drawer-metric-val">Identical amino acid change known</span></div>' : '') +
          (v.clinvarId ? '<div class="drawer-metric-row"><span class="drawer-metric-lbl">ClinVar Variation:</span><span class="drawer-metric-val"><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/' + v.clinvarId + '/" target="_blank" rel="noopener">VCV#' + v.clinvarId + ' ↗</a></span></div>' : '') +
        '</div>' +

        // Box 4: Research & GWAS Links
        '<div class="variant-drawer-box">' +
          '<h4>Research Studies</h4>' +
          (v.studies && v.studies.length
            ? v.studies.map((s) =>
                '<div style="margin-bottom:4px;">' +
                  '<div style="font-weight:700;font-size:11.5px;color:var(--teal-dark);">' + s.title + '</div>' +
                  '<div style="font-size:11px;color:var(--slate);">' + s.description + '</div>' +
                  '<div style="font-size:10.5px;color:var(--slate);margin-top:2px;">' + s.source + ' · ' + s.genotypeRelevance + '</div>' +
                '</div>'
              ).join("")
            : '<div style="font-size:11.5px;color:var(--slate);">No individual GWAS study correlations published for this exact SNP locus.</div>') +
        '</div>' +
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
            '<span class="study-card__variant" data-gene="' + (v.gene || '') + '" data-vid="' + v.id + '">' + s.title + '</span>' +
            '<span class="study-card__condition">' + s.condition + '</span>' +
          '</div>' +
          '<div class="study-card__finding" style="font-weight:600;margin-top:3px;">' + s.finding + '</div>' +
          '<div style="font-size:12px;color:var(--slate);margin:4px 0;">' + s.description + '</div>' +
          '<div class="study-card__foot">' +
            '<span>' + s.source + '</span>' +
            '<span>' + s.genotypeRelevance + (s.pValue && s.pValue !== 'N/A' ? ' · p=' + s.pValue : '') + '</span>' +
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

    const width = 720, height = 180, padX = 40;
    const svgns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgns, "svg");
    svg.setAttribute("class", "genome-browser-svg");
    svg.setAttribute("width", width);
    svg.setAttribute("height", height);
    svg.setAttribute("viewBox", "0 0 " + width + " " + height);

    const axisY = 120;
    const axisLine = document.createElementNS(svgns, "line");
    axisLine.setAttribute("x1", padX); axisLine.setAttribute("y1", axisY);
    axisLine.setAttribute("x2", width - padX); axisLine.setAttribute("y2", axisY);
    axisLine.setAttribute("stroke", "var(--line-strong)");
    axisLine.setAttribute("stroke-width", "2");
    svg.appendChild(axisLine);

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

    const riskGrid = byId("organ-risk-grid");
    if (riskGrid && typeof ORGAN_RISK_MATRIX !== "undefined") {
      riskGrid.innerHTML = ORGAN_RISK_MATRIX.map((m) =>
        '<div class="risk-matrix-card">' +
          '<div class="risk-matrix-card__head">' +
            '<div class="risk-matrix-card__title"><span>' + m.icon + '</span><span>' + m.system + '</span></div>' +
            '<span class="risk-matrix-card__tier ' + m.riskTier + '">' + m.riskTier + ' RISK</span>' +
          '</div>' +
          '<div class="risk-matrix-card__body">' +
            '<div><strong>Primary Pathway:</strong> ' + m.pathway + '</div>' +
            '<div style="margin-top:3px;"><strong>Polygenic Percentile:</strong> ' + m.prsPercentile + 'th percentile</div>' +
            (m.concernGenes && m.concernGenes.length
              ? '<div class="risk-matrix-card__genes">' +
                m.concernGenes.map((g) => '<span class="risk-matrix-card__gene">' + g + '</span>').join("") +
                '</div>'
              : '') +
          '</div>' +
        '</div>'
      ).join("");
    }

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
      '<div class="pub-card__authors">' + (p.authors || "Consortium") + ' · <span class="mono">' + (p.journal || "PubMed") + ' (' + (p.year || "2023") + ')</span></div>' +
      (p.relevance ? '<div class="pub-card__relevance" style="font-size:12px;color:var(--ink);margin:4px 0;">' + p.relevance + '</div>' : '') +
      '<a class="pub-card__link" href="' + (p.url || 'https://pubmed.ncbi.nlm.nih.gov/' + p.pmid + '/') + '" target="_blank" rel="noopener">PMID ' + p.pmid + ' ↗</a>' +
      '</div>'
    );
  }

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

    ["filter-gene", "filter-clinvar", "filter-phase", "filter-zygosity"].forEach((id) => {
      const el = byId(id);
      if (el) el.addEventListener("change", renderVariantsTable);
    });

    const geneSearch = byId("genes-search");
    if (geneSearch) geneSearch.addEventListener("input", renderGenesTable);

    // Initial render
    renderLeftPanel();

    const filtered = filteredOntology();
    if (filtered.length > 0) {
      selectCategory(filtered[0]);
    }
  });

})();
