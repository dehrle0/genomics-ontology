/**
 * GENOMIC COMMAND CENTER & VARIANT ATLAS — APPLICATION ENGINE
 * Pure SVG responsive charts, Screen-1 Executive Wisdom Ribbon, Multi-Engine Concordance
 * Quadrant, Intra-Gene Cis/Trans Phase Auditor, Collapsible A–Z Gene Index,
 * AlphaGenome 18-Modality Sparkbars, and Interactive SVG Gene Relationship Network Modal.
 */

(function () {
  "use strict";

  // Application State
  const state = {
    currentPatient: "Daniel_Ehrle",
    currentReport: "visual_explorer",
    geneScope: "all", // "all", "tier1_top15", "tier1_only", "top10", "top15", "top25"
    searchTerm: "",
    categoryFilter: null,
    quadrantFilter: null,
    phaseFilter: null,
    barScale: "linear", // "linear" or "log10"
    sidebarCollapsed: false,
    sidebarSort: "alpha", // "alpha" (A-Z) or "score" (Impact)
    sidebarSearch: "",
    sidebarLetter: "ALL",
    selectedGeneSymbol: null,
    openGenes: new Set(),
    openVariants: new Set(),
    cachedNcbiSummaries: {}
  };

  // Helper: DOM Selector
  const $ = (id) => document.getElementById(id);

  // Initialize Application
  window.addEventListener("DOMContentLoaded", () => {
    initCohortSwitcher();
    initReportSelector();
    initSearchFilter();
    initToolbarActions();
    initSidebarControls();
    initNetworkModal();
    loadPatientData(state.currentPatient);
  });

  // ---------------------------------------------------------------------------
  // 1. COHORT SWITCHER & HEADER SETUP
  // ---------------------------------------------------------------------------
  function initCohortSwitcher() {
    const btns = document.querySelectorAll(".cohort-btn");
    btns.forEach((btn) => {
      btn.addEventListener("click", () => {
        const patientId = btn.getAttribute("data-patient");
        if (patientId === state.currentPatient) return;

        btns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");

        state.currentPatient = patientId;
        state.openGenes.clear();
        state.openVariants.clear();
        state.categoryFilter = null;
        state.quadrantFilter = null;
        state.phaseFilter = null;
        state.selectedGeneSymbol = null;
        loadPatientData(patientId);
      });
    });
  }

  function loadPatientData(patientId) {
    if (!window.COHORT_DATA) {
      console.error("window.COHORT_DATA is missing.");
      return;
    }
    // Fallback if standalone single-patient HTML doesn't have the clicked cohort
    if (!window.COHORT_DATA[patientId]) {
      const available = Object.keys(window.COHORT_DATA);
      if (available.length > 0) {
        patientId = available[0];
        state.currentPatient = patientId;
      } else {
        return;
      }
    }
    const pData = window.COHORT_DATA[patientId];

    // Update Header Meta & Audience
    const metaEl = $("topbar-meta");
    if (metaEl) {
      const dpInfo = pData.meta?.vcfDepthTelemetry?.mean_depth_wgs || "Phased WGS";
      metaEl.textContent = `${pData.patientName} · ${pData.meta.totalVariants} variants · ${pData.meta.totalGenes} genes · Mean DP: ${dpInfo}`;
    }
    const audEl = $("topbar-audience");
    if (audEl && pData.meta?.targetAudience) {
      audEl.textContent = `Target Audience: ${pData.meta.targetAudience}`;
    }

    // Render Screen-1 Executive Wisdom Ribbon & Concordance Matrix
    renderWisdomSection(pData);

    // Render 3 Rows of Landscape Charts (Ontology, Organ Systems, Pharmacogenomics)
    renderAllLandscapeCharts(pData.charts);

    // Render Current Report View + Left Gene Sidebar
    renderActiveReport();
  }

  // ---------------------------------------------------------------------------
  // 2. SCREEN-1 EXECUTIVE WISDOM RIBBON & CONCORDANCE / PHASE MATRIX
  // ---------------------------------------------------------------------------
  function renderWisdomSection(pData) {
    const ribbonMount = $("wisdom-ribbon-mount");
    const matrixMount = $("wisdom-matrix-mount");
    if (!ribbonMount || !matrixMount) return;

    const wisdom = pData.wisdom || {};
    const imp = wisdom.imperatives || {};
    const cisTransList = wisdom.cisTransAudit || [];
    const quads = wisdom.concordanceQuadrants || {};
    const depth = pData.meta?.vcfDepthTelemetry || {};

    const rxList = imp.rxHighlights || imp.rxGuardrails || [];
    const monoList = imp.monogenicHighlights || imp.monogenicfindings || [];
    const protList = imp.protectiveHighlights || imp.protectiveOffsets || [];
    const transGenes = cisTransList.filter((g) => (g.phaseState || "").includes("Trans"));
    const cisGenes = cisTransList.filter((g) => (g.phaseState || "").includes("Cis"));
    const transCount = imp.transCompoundCandidatesCount !== undefined ? imp.transCompoundCandidatesCount : transGenes.length;
    const cisCount = imp.cisHaplotypesCount !== undefined ? imp.cisHaplotypesCount : cisGenes.length;

    // 4-Card Executive Wisdom Ribbon
    ribbonMount.innerHTML = `
      <div class="wisdom-ribbon-grid">
        <!-- Card 1: Critical Rx & Perioperative Guardrails -->
        <div class="wisdom-card wisdom-card--rx" id="wisdom-card-rx" title="Click to view Pharmacogenomic &amp; Tier-1 Prescribing Guardrails">
          <div class="wisdom-card__top">
            <span class="wisdom-card__label">🚨 Critical Rx &amp; Perioperative Guardrails</span>
            <span class="wisdom-card__kpi" style="color:var(--concern);">${imp.rxGuardrailsCount !== undefined ? imp.rxGuardrailsCount : rxList.length}</span>
          </div>
          <div class="wisdom-card__body">
            CPIC/DPWG actionable drug-gene metabolizer alerts (Avoid, Reduce Dose, Increase Dose, Caution) requiring dosage adjustment or alternative selection.
          </div>
          <div class="wisdom-card__pills">
            ${rxList.slice(0, 4).map((r) => `<span class="wisdom-pill" style="border-color:#fecaca;background:#fff1f2;color:#b91c1c;">${r.gene}: ${r.drug} (${r.action || "PGx"})</span>`).join("") || `<span class="wisdom-pill">Standard PGx Panel</span>`}
          </div>
        </div>

        <!-- Card 2: Monogenic & High-Penetrance Findings -->
        <div class="wisdom-card wisdom-card--monogenic" id="wisdom-card-mono" title="Click to filter to ClinVar Pathogenic / Likely Pathogenic findings">
          <div class="wisdom-card__top">
            <span class="wisdom-card__label">🧬 Monogenic &amp; Carrier Findings</span>
            <span class="wisdom-card__kpi" style="color:#b45309;">${imp.monogenicCount !== undefined ? imp.monogenicCount : monoList.length}</span>
          </div>
          <div class="wisdom-card__body">
            ClinVar Pathogenic and Likely Pathogenic alleles curated across Mendelian &amp; high-penetrance clinical loci.
          </div>
          <div class="wisdom-card__pills">
            ${monoList.slice(0, 5).map((m) => `<span class="wisdom-pill" style="border-color:#fde68a;background:#fffbeb;color:#92400e;" data-jump-gene="${m.gene}">${m.gene} (${m.zygosity || "Het"})</span>`).join("") || `<span class="wisdom-pill">View ClinVar Catalog</span>`}
          </div>
        </div>

        <!-- Card 3: Intra-Gene Cis/Trans Phase Auditor -->
        <div class="wisdom-card wisdom-card--phase" id="wisdom-card-phase" title="Click to filter to Multi-Hit Trans vs. Cis Haplotype Genes">
          <div class="wisdom-card__top">
            <span class="wisdom-card__label">⚖️ Intra-Gene Cis / Trans Phase Auditor</span>
            <span class="wisdom-card__kpi" style="color:#1d4ed8;">${transCount} Trans · ${cisCount} Cis</span>
          </div>
          <div class="wisdom-card__body">
            Parental &amp; read-backed phasing distinguishing multi-variant <strong>Trans</strong> genes from uni-allelic <strong>Cis</strong> haplotypes (1 intact WT copy).
          </div>
          <div class="wisdom-card__pills">
            ${transGenes.slice(0, 3).map((g) => `<span class="wisdom-pill" style="border-color:#fecaca;background:#fef2f2;color:#b91c1c;" data-jump-gene="${g.gene}">${g.gene} (Trans)</span>`).join("")}
            ${cisGenes.slice(0, 2).map((g) => `<span class="wisdom-pill" style="border-color:#bfdbfe;background:#eff6ff;color:#1d4ed8;" data-jump-gene="${g.gene}">${g.gene} (Cis)</span>`).join("")}
          </div>
        </div>

        <!-- Card 4: Protective & Longevity Offsets -->
        <div class="wisdom-card wisdom-card--protect" id="wisdom-card-protect" title="Click to filter to Protective &amp; Resilience Alleles">
          <div class="wisdom-card__top">
            <span class="wisdom-card__label">🛡️ Protective &amp; Physiological Offsets</span>
            <span class="wisdom-card__kpi" style="color:var(--protect);">${imp.protectiveCount !== undefined ? imp.protectiveCount : protList.length}</span>
          </div>
          <div class="wisdom-card__body">
            Curated protective, resilience, and risk-attenuating alleles buffering polygenic and metabolic pathway susceptibility.
          </div>
          <div class="wisdom-card__pills">
            ${protList.slice(0, 5).map((p) => `<span class="wisdom-pill" style="border-color:#bbf7d0;background:#f0fdf4;color:#166534;" data-jump-gene="${p.gene}">${p.gene} (${p.variant || p.id || "Protective"})</span>`).join("") || `<span class="wisdom-pill">View Protective Loci</span>`}
          </div>
        </div>
      </div>
    `;

    // Matrix Panel: Sequencing Depth + Multi-Engine Concordance Quadrant + Cis/Trans Phase Auditor
    const qOrder = [
      "q1_concordant_actionable",
      "q2_ai_rescue_vus",
      "q3_clinical_mod_low_ai",
      "q4_protective_offset"
    ];

    matrixMount.innerHTML = `
      <div class="wisdom-matrix-card">
        <div class="wisdom-matrix__header">
          <div class="wisdom-matrix__title">
            <span>🧭 Epistemic Wisdom Synthesis: Multi-Engine Concordance &amp; Sequencing Depth</span>
            ${state.quadrantFilter || state.phaseFilter ? `<button class="btn-secondary" id="btn-clear-wisdom-filter" style="padding:2px 8px;font-size:11px;color:var(--concern);border-color:#fecaca;">✕ Clear Active Wisdom Filter</button>` : ""}
          </div>
          <div class="depth-telemetry-strip">
            <span class="depth-pill" title="Empirical mean coverage across called PASS VCF loci">Mean WGS Depth: <strong>${depth.mean_depth_called_loci || depth.mean_depth_wgs || "50.1x"}</strong> (Median <strong>${depth.median_depth_called_loci || "50.0x"}</strong>)</span>
            <span class="depth-pill" title="Reference genome assembly">Reference: <strong>GRCh38.p14 Phased WGS</strong></span>
          </div>
        </div>

        <div class="wisdom-matrix-grid">
          <!-- Left: Multi-Engine Concordance Quadrant -->
          <div>
            <div style="font-size:11.5px;font-weight:700;color:var(--slate);margin-bottom:6px;display:flex;justify-content:space-between;">
              <span>MULTI-ENGINE CONCORDANCE QUADRANT (CLINVAR vs. AI FUNCTIONAL IMPACT)</span>
              <span style="font-weight:500;font-size:10.5px;">Click quadrant to filter explorer</span>
            </div>
            <div class="concordance-grid">
              ${qOrder.map((qKey) => {
                const qObj = quads[qKey] || { label: qKey, count: 0, desc: "" };
                const isActive = state.quadrantFilter === qKey;
                return `
                  <div class="concordance-tile ${isActive ? "is-active" : ""}" data-quadrant="${qKey}">
                    <div class="concordance-tile__top">
                      <span class="concordance-tile__name">${qObj.label}</span>
                      <span class="concordance-tile__count">${qObj.count}</span>
                    </div>
                    <div class="concordance-tile__desc">${qObj.desc}</div>
                  </div>
                `;
              }).join("")}
            </div>
          </div>

          <!-- Right: Intra-Gene Cis/Trans Phase Auditor -->
          <div class="phase-audit-box">
            <div style="display:flex;justify-content:space-between;align-items:center;">
              <span style="font-size:11.5px;font-weight:700;color:var(--slate);">INTRA-GENE CIS / TRANS PHASE AUDITOR (MULTI-HIT GENES)</span>
              <div style="display:flex;gap:4px;">
                <button class="btn-secondary" data-phase-filter="Trans" style="padding:1px 7px;font-size:10px;${state.phaseFilter === "Trans" ? "background:#fef2f2;border-color:#ef4444;color:#b91c1c;" : ""}">Filter Trans (${transCount})</button>
                <button class="btn-secondary" data-phase-filter="Cis" style="padding:1px 7px;font-size:10px;${state.phaseFilter === "Cis" ? "background:#eff6ff;border-color:#3b82f6;color:#1d4ed8;" : ""}">Filter Cis (${cisCount})</button>
              </div>
            </div>
            <div style="font-size:11px;color:#475569;">
              Click any gene chip below to jump directly to its phased allelic table (<span style="color:#b91c1c;font-weight:700;">Red = Trans Multi-Hit</span> · <span style="color:#1d4ed8;font-weight:700;">Blue = Cis Haplotype</span>):
            </div>
            <div class="phase-audit-list">
              ${cisTransList.slice(0, 36).map((ag) => {
                const isTrans = (ag.phaseState || "").includes("Trans");
                return `<span class="phase-gene-chip ${isTrans ? "phase-gene-chip--trans" : "phase-gene-chip--cis"}" data-jump-gene="${ag.gene}" title="${ag.gene}: ${ag.phaseState} (${ag.totalVariants} variants: ${ag.maternalCount} Maternal, ${ag.paternalCount} Paternal) — ${ag.clinicalNote || ""}">${ag.gene} (${isTrans ? "Trans" : "Cis"} · ${ag.totalVariants})</span>`;
              }).join("") || `<span style="font-size:11px;color:var(--slate);">Single-variant locus distribution.</span>`}
            </div>
          </div>
        </div>
      </div>
    `;

    // Attach Wisdom Event Listeners
    const rxCard = $("wisdom-card-rx");
    if (rxCard) {
      rxCard.addEventListener("click", (e) => {
        if (e.target.hasAttribute("data-jump-gene")) return;
        const sel = $("report-selector");
        if (sel) {
          sel.value = "ai_clinical";
          state.currentReport = "ai_clinical";
        }
        renderActiveReport();
      });
    }

    const monoCard = $("wisdom-card-mono");
    if (monoCard) {
      monoCard.addEventListener("click", (e) => {
        if (e.target.hasAttribute("data-jump-gene")) return;
        const sel = $("report-selector");
        if (sel) {
          sel.value = "clinvar";
          state.currentReport = "clinvar";
        }
        state.quadrantFilter = null;
        renderWisdomSection(pData);
        renderActiveReport();
      });
    }

    const phaseCard = $("wisdom-card-phase");
    if (phaseCard) {
      phaseCard.addEventListener("click", (e) => {
        if (e.target.hasAttribute("data-jump-gene")) return;
        state.phaseFilter = state.phaseFilter === "Trans" ? null : "Trans";
        renderWisdomSection(pData);
        renderActiveReport();
      });
    }

    const protCard = $("wisdom-card-protect");
    if (protCard) {
      protCard.addEventListener("click", (e) => {
        if (e.target.hasAttribute("data-jump-gene")) return;
        const sel = $("report-selector");
        if (sel) {
          sel.value = "clinvar";
          state.currentReport = "clinvar";
        }
        state.searchTerm = "protective";
        const sInput = $("search-input");
        if (sInput) sInput.value = "protective";
        renderActiveReport();
      });
    }

    // Quadrant tile clicks
    matrixMount.querySelectorAll("[data-quadrant]").forEach((tile) => {
      tile.addEventListener("click", () => {
        const qk = tile.getAttribute("data-quadrant");
        state.quadrantFilter = state.quadrantFilter === qk ? null : qk;
        renderWisdomSection(pData);
        renderActiveReport();
      });
    });

    // Phase filter buttons
    matrixMount.querySelectorAll("[data-phase-filter]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const pf = btn.getAttribute("data-phase-filter");
        state.phaseFilter = state.phaseFilter === pf ? null : pf;
        renderWisdomSection(pData);
        renderActiveReport();
      });
    });

    const clearWisBtn = $("btn-clear-wisdom-filter");
    if (clearWisBtn) {
      clearWisBtn.addEventListener("click", () => {
        state.quadrantFilter = null;
        state.phaseFilter = null;
        renderWisdomSection(pData);
        renderActiveReport();
      });
    }

    // Direct jump-to-gene pills
    document.querySelectorAll("[data-jump-gene]").forEach((pill) => {
      pill.addEventListener("click", (e) => {
        e.stopPropagation();
        const sym = pill.getAttribute("data-jump-gene");
        jumpToGeneCard(sym);
      });
    });
  }

  function jumpToGeneCard(symbol) {
    if (!symbol) return;
    if (state.currentReport === "ai_deep" || state.currentReport === "ai_clinical") {
      state.currentReport = "visual_explorer";
      const sel = $("report-selector");
      if (sel) sel.value = "visual_explorer";
    }
    state.categoryFilter = null;
    state.quadrantFilter = null;
    state.phaseFilter = null;
    state.searchTerm = "";
    const sInput = $("search-input");
    if (sInput) sInput.value = "";

    state.openGenes.add(symbol);
    state.selectedGeneSymbol = symbol;
    renderActiveReport();

    setTimeout(() => {
      const cardEl = $(`gene-card-${symbol}`);
      if (cardEl) {
        cardEl.scrollIntoView({ behavior: "smooth", block: "center" });
      }
    }, 60);
  }

  // ---------------------------------------------------------------------------
  // 3. SVG CHART ENGINE (PANEL A: STACKED BAR | PANELS B & C: PIE WITH CALLOUTS)
  // ---------------------------------------------------------------------------
  function renderAllLandscapeCharts(charts) {
    if (!charts) return;

    // Row 1: Groupings by Ontology
    if (charts.row1_ontology) {
      renderStackedBar("chart-row1-bar", charts.row1_ontology.bar, "Diagnostic vs Secondary by HPO Domain", "Variants", true);
      renderPieChart("chart-row1-pie-prior", charts.row1_ontology.piePrior, "Presentation prior to RGS");
      renderPieChart("chart-row1-pie-post", charts.row1_ontology.piePost, "Diagnosis Category after RGS");
    }

    // Row 2: Groupings by Organ System (Restored)
    if ($("chart-row2-bar") && charts.row2_organ) {
      renderStackedBar("chart-row2-bar", charts.row2_organ.bar, "Actionability by Anatomical Organ System", "Variants", true);
      renderPieChart("chart-row2-pie-prior", charts.row2_organ.piePrior, "Organ Signals prior to Resolution");
      renderPieChart("chart-row2-pie-post", charts.row2_organ.piePost, "Primary Confirmed Organ Attribution");
    }

    // Row 3: Groupings by Pharmacogenomics
    if (charts.row3_pharma) {
      renderStackedBar("chart-row3-bar", charts.row3_pharma.bar, "Actionability by Therapeutic Drug Class", "Drug Cards", false);
      renderPieChart("chart-row3-pie-prior", charts.row3_pharma.piePrior, "Metabolizer Phenotypes");
      renderPieChart("chart-row3-pie-post", charts.row3_pharma.piePost, "Clinical Action Categories");
    }
  }

  function renderStackedBar(containerId, data, title, yUnit, allowScaleToggle = false) {
    const mount = $(containerId);
    if (!mount || !data) return;
    mount.innerHTML = "";

    const isLog = state.barScale === "log10";
    const W = 520;
    const H = 270;
    const top = 40;
    const bottom = 72;
    const left = 50;
    const right = 20;
    const chartW = W - left - right;
    const chartH = H - top - bottom;

    const maxRaw = Math.max(...data.map((d) => d.total), 1);
    const toScale = (v) => (isLog ? (v > 0 ? Math.log10(v + 1) : 0) : v);
    const maxScaled = toScale(maxRaw) * 1.15;

    const n = Math.max(data.length, 1);
    const colW = chartW / n;
    const barW = Math.max(colW * 0.58, 14);

    let html = "";
    if (allowScaleToggle) {
      html += `
        <div style="display:flex;justify-content:flex-end;width:100%;margin-bottom:4px;">
          <div class="scale-toggle">
            <button class="scale-btn ${!isLog ? "active" : ""}" data-scale-target="linear">Linear</button>
            <button class="scale-btn ${isLog ? "active" : ""}" data-scale-target="log10">Log₁₀</button>
          </div>
        </div>
      `;
    }

    let svg = `<svg viewBox="0 0 ${W} ${H}" class="chart-svg" xmlns="http://www.w3.org/2000/svg">`;

    if (!isLog) {
      const numTicks = 5;
      for (let i = 0; i <= numTicks; i++) {
        const val = Math.round((maxRaw / numTicks) * i);
        const y = top + chartH - (val / (maxRaw * 1.15)) * chartH;
        svg += `<line x1="${left}" y1="${y}" x2="${left + chartW}" y2="${y}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3" />`;
        svg += `<text x="${left - 8}" y="${y + 3.5}" text-anchor="end" font-size="10" font-family="system-ui" fill="#64748b">${val}</text>`;
      }
    } else {
      const logTicks = [1, 10, 50, 100, 250, 500, 1000].filter((t) => t <= maxRaw * 1.2);
      logTicks.forEach((t) => {
        const y = top + chartH - (toScale(t) / maxScaled) * chartH;
        svg += `<line x1="${left}" y1="${y}" x2="${left + chartW}" y2="${y}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3" />`;
        svg += `<text x="${left - 8}" y="${y + 3.5}" text-anchor="end" font-size="9.5" font-family="system-ui" fill="#64748b">${t}</text>`;
      });
    }

    svg += `<line x1="${left}" y1="${top}" x2="${left}" y2="${top + chartH}" stroke="#17222b" stroke-width="1.5" />`;
    svg += `<text transform="translate(14, ${top + chartH / 2}) rotate(-90)" text-anchor="middle" font-size="10.5" font-weight="600" fill="#475569">${yUnit} ${isLog ? "(Log₁₀)" : ""}</text>`;
    svg += `<line x1="${left}" y1="${top + chartH}" x2="${left + chartW}" y2="${top + chartH}" stroke="#17222b" stroke-width="1.5" />`;

    data.forEach((d, i) => {
      const barX = left + i * colW + (colW - barW) / 2;
      const totalH = (toScale(d.total) / maxScaled) * chartH;
      const actH = (toScale(d.actionable) / maxScaled) * chartH;
      const secH = Math.max(totalH - actH, 0);

      const actY = top + chartH - totalH;
      const secY = top + chartH - secH;

      svg += `<text x="${barX + barW / 2}" y="${actY - 6}" text-anchor="middle" font-size="10" font-weight="700" fill="#17222b">${d.fraction}</text>`;

      if (actH > 0) {
        svg += `<rect x="${barX}" y="${actY}" width="${barW}" height="${actH}" fill="#B3261E" rx="1" style="cursor:pointer;" data-category="${d.category}">
          <title>${d.category}: ${d.actionable} Actionable / ${d.total} Total (${d.ratePct}%)</title>
        </rect>`;
      }

      if (secH > 0) {
        svg += `<rect x="${barX}" y="${secY}" width="${barW}" height="${secH}" fill="#17222B" rx="1" style="cursor:pointer;" data-category="${d.category}">
          <title>${d.category}: ${d.secondary} Secondary / ${d.total} Total</title>
        </rect>`;
      }

      const lblX = barX + barW / 2;
      const lblY = top + chartH + 12;
      const shortCat = d.category.length > 14 ? d.category.slice(0, 13) + "…" : d.category;
      svg += `<text transform="translate(${lblX}, ${lblY}) rotate(-35)" text-anchor="end" font-size="10" font-weight="600" fill="#334155" style="cursor:pointer;" data-category="${d.category}">${shortCat}</text>`;
    });

    const legY = H - 8;
    svg += `
      <g transform="translate(${W / 2 - 110}, ${legY})">
        <rect x="0" y="-8" width="10" height="10" fill="#17222B" rx="1" />
        <text x="14" y="0" font-size="10" font-weight="600" fill="#475569">Secondary / VUS</text>
        <rect x="120" y="-8" width="10" height="10" fill="#B3261E" rx="1" />
        <text x="134" y="0" font-size="10" font-weight="600" fill="#475569">Actionable / Diagnostic</text>
      </g>
    `;

    svg += `</svg>`;
    mount.innerHTML = html + svg;

    mount.querySelectorAll("[data-category]").forEach((el) => {
      el.addEventListener("click", () => {
        const cat = el.getAttribute("data-category");
        setCategoryFilter(cat);
      });
    });

    mount.querySelectorAll("[data-scale-target]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const sc = btn.getAttribute("data-scale-target");
        if (state.barScale !== sc) {
          state.barScale = sc;
          const pData = window.COHORT_DATA[state.currentPatient];
          if (pData) renderAllLandscapeCharts(pData.charts);
        }
      });
    });
  }

  function renderPieChart(containerId, slices, title) {
    const mount = $(containerId);
    if (!mount || !slices) return;
    mount.innerHTML = "";

    const W = 380;
    const H = 265;
    const cx = 190;
    const cy = 132;
    const radius = 72;

    const total = slices.reduce((acc, s) => acc + s.count, 0) || 1;
    let svg = `<svg viewBox="0 0 ${W} ${H}" class="chart-svg" xmlns="http://www.w3.org/2000/svg">`;

    let currentAngle = -Math.PI / 2;
    const rawCallouts = [];
    slices.forEach((s) => {
      if (s.count <= 0) return;
      const sliceAngle = (s.count / total) * 2 * Math.PI;
      const startAngle = currentAngle;
      const endAngle = currentAngle + sliceAngle;
      currentAngle = endAngle;

      const x1 = cx + radius * Math.cos(startAngle);
      const y1 = cy + radius * Math.sin(startAngle);
      const x2 = cx + radius * Math.cos(endAngle);
      const y2 = cy + radius * Math.sin(endAngle);

      const largeArc = sliceAngle > Math.PI ? 1 : 0;
      const pathD = sliceAngle >= 2 * Math.PI - 0.001
        ? `M ${cx} ${cy - radius} A ${radius} ${radius} 0 1 1 ${cx - 0.01} ${cy - radius} Z`
        : `M ${cx} ${cy} L ${x1} ${y1} A ${radius} ${radius} 0 ${largeArc} 1 ${x2} ${y2} Z`;

      svg += `<path d="${pathD}" fill="${s.color}" stroke="#ffffff" stroke-width="1.2" style="cursor:pointer;" data-label="${s.label}">
        <title>${s.label}: ${s.pct}% (${s.count}/${total})</title>
      </path>`;

      const midAngle = startAngle + sliceAngle / 2;
      rawCallouts.push({
        slice: s,
        midAngle: midAngle,
        innerX: cx + radius * Math.cos(midAngle),
        innerY: cy + radius * Math.sin(midAngle),
        outerX: cx + (radius + 16) * Math.cos(midAngle),
        outerY: cy + (radius + 16) * Math.sin(midAngle)
      });
    });

    const rightSide = rawCallouts.filter((c) => Math.cos(c.midAngle) >= 0);
    const leftSide = rawCallouts.filter((c) => Math.cos(c.midAngle) < 0);

    rightSide.sort((a, b) => a.outerY - b.outerY);
    leftSide.sort((a, b) => a.outerY - b.outerY);

    function relaxVerticalPositions(list, boxH, gap, minY, maxY) {
      if (list.length === 0) return;
      const minStep = boxH + gap;
      for (let i = 1; i < list.length; i++) {
        if (list[i].outerY < list[i - 1].outerY + minStep) {
          list[i].outerY = list[i - 1].outerY + minStep;
        }
      }
      if (list[list.length - 1].outerY > maxY) {
        list[list.length - 1].outerY = maxY;
        for (let i = list.length - 2; i >= 0; i--) {
          if (list[i].outerY > list[i + 1].outerY - minStep) {
            list[i].outerY = list[i + 1].outerY - minStep;
          }
        }
      }
      if (list[0].outerY < minY) {
        list[0].outerY = minY;
        for (let i = 1; i < list.length; i++) {
          if (list[i].outerY < list[i - 1].outerY + minStep) {
            list[i].outerY = list[i - 1].outerY + minStep;
          }
        }
      }
    }

    const boxW = 88;
    const boxH = 24;
    relaxVerticalPositions(rightSide, boxH, 4, 18, H - 22);
    relaxVerticalPositions(leftSide, boxH, 4, 18, H - 22);

    const adjustedCallouts = [...rightSide, ...leftSide];
    adjustedCallouts.forEach((c) => {
      const isRight = Math.cos(c.midAngle) >= 0;
      const elbowX = isRight ? c.outerX + 14 : c.outerX - 14;
      const boxX = isRight ? elbowX + 2 : elbowX - boxW - 2;
      const boxY = c.outerY - boxH / 2;

      svg += `<polyline points="${c.innerX},${c.innerY} ${c.outerX},${c.outerY} ${elbowX},${c.outerY}" stroke="#64748b" stroke-width="1" fill="none" />`;
      svg += `<g style="cursor:pointer;" data-label="${c.slice.label}">`;
      svg += `<rect x="${boxX}" y="${boxY}" width="${boxW}" height="${boxH}" fill="#ffffff" stroke="#94a3b8" stroke-width="0.8" rx="3" filter="drop-shadow(0 1px 2px rgba(0,0,0,0.06))" />`;

      const lblText = c.slice.label.length > 13 ? c.slice.label.slice(0, 12) + "…" : c.slice.label;
      svg += `<text x="${boxX + boxW / 2}" y="${boxY + 10}" text-anchor="middle" font-size="8.5" font-weight="700" fill="#1e293b">${lblText}</text>`;
      svg += `<text x="${boxX + boxW / 2}" y="${boxY + 20}" text-anchor="middle" font-size="8.5" font-weight="600" fill="#0f766e">${c.slice.pct}%, <tspan fill="#64748b">${c.slice.count}</tspan></text>`;
      svg += `</g>`;
    });

    svg += `</svg>`;

    let legendHtml = `<div class="pie-slices-legend">`;
    slices.forEach((s) => {
      if (s.count <= 0) return;
      legendHtml += `
        <span class="pie-slice-pill" style="cursor:pointer;" data-label="${s.label}">
          <span class="pie-slice-dot" style="background:${s.color};"></span>
          <span class="pie-slice-name">${s.label}:</span>
          <strong class="pie-slice-pct">${s.pct}%</strong>
          <span class="pie-slice-fraction mono">(${s.count})</span>
        </span>
      `;
    });
    legendHtml += `</div>`;

    mount.innerHTML = svg + legendHtml;

    mount.querySelectorAll("[data-label]").forEach((el) => {
      el.addEventListener("click", () => {
        const lbl = el.getAttribute("data-label");
        setCategoryFilter(lbl);
      });
    });
  }

  function setCategoryFilter(category) {
    state.categoryFilter = state.categoryFilter === category ? null : category;
    renderActiveReport();
  }

  // ---------------------------------------------------------------------------
  // 4. COLLAPSIBLE LEFT COLUMN (A-Z & IMPACT GENE NAVIGATOR)
  // ---------------------------------------------------------------------------
  function initSidebarControls() {
    const toggleBtn = $("btn-toggle-sidebar");
    const ws = $("explorer-workspace");
    if (toggleBtn && ws) {
      toggleBtn.addEventListener("click", () => {
        state.sidebarCollapsed = !state.sidebarCollapsed;
        ws.classList.toggle("sidebar-collapsed", state.sidebarCollapsed);
      });
    }

    const sortAlpha = $("sort-alpha");
    const sortScore = $("sort-score");
    if (sortAlpha && sortScore) {
      sortAlpha.addEventListener("click", () => {
        state.sidebarSort = "alpha";
        sortAlpha.classList.add("active");
        sortScore.classList.remove("active");
        renderActiveReport();
      });
      sortScore.addEventListener("click", () => {
        state.sidebarSort = "score";
        sortScore.classList.add("active");
        sortAlpha.classList.remove("active");
        renderActiveReport();
      });
    }

    const sideSearch = $("sidebar-gene-search");
    if (sideSearch) {
      sideSearch.addEventListener("input", (e) => {
        state.sidebarSearch = e.target.value.trim().toLowerCase();
        renderActiveReport();
      });
    }
  }

  function renderSidebarGeneNavigator(processedGenes) {
    const countEl = $("sidebar-gene-count");
    const azMount = $("sidebar-az-bar");
    const listMount = $("sidebar-gene-list");
    if (!listMount) return;

    // Build A-Z bar
    const letters = ["ALL", ..."ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("")];
    if (azMount) {
      azMount.innerHTML = letters.map((L) => `
        <button class="az-btn ${state.sidebarLetter === L ? "active" : ""}" data-letter="${L}" style="${L === "ALL" ? "width:auto;padding:0 5px;" : ""}">${L}</button>
      `).join("");

      azMount.querySelectorAll(".az-btn").forEach((b) => {
        b.addEventListener("click", () => {
          state.sidebarLetter = b.getAttribute("data-letter");
          renderSidebarGeneNavigator(processedGenes);
        });
      });
    }

    // Clone and sort for sidebar
    let navGenes = [...processedGenes];
    if (state.sidebarLetter && state.sidebarLetter !== "ALL") {
      navGenes = navGenes.filter((g) => (g.symbol || "").toUpperCase().startsWith(state.sidebarLetter));
    }
    if (state.sidebarSearch) {
      navGenes = navGenes.filter((g) =>
        (g.symbol || "").toLowerCase().includes(state.sidebarSearch) ||
        (g.name || "").toLowerCase().includes(state.sidebarSearch)
      );
    }

    if (state.sidebarSort === "alpha") {
      navGenes.sort((a, b) => (a.symbol || "").localeCompare(b.symbol || ""));
    } else {
      navGenes.sort((a, b) => b.geneScore - a.geneScore);
    }

    if (countEl) countEl.textContent = navGenes.length;

    if (navGenes.length === 0) {
      listMount.innerHTML = `<div style="padding:12px;font-size:11px;color:var(--slate);text-align:center;">No matching genes in index.</div>`;
      return;
    }

    listMount.innerHTML = navGenes.map((g) => {
      const hasPath = g.variants.some((v) => v.tier === "Tier1" || (v.clinvar || "").toLowerCase().includes("pathogenic"));
      const dotColor = hasPath ? "#b3261e" : (g.phaseAudit?.phaseStatus?.includes("Trans") ? "#d97706" : "#0f766e");
      const isSel = state.selectedGeneSymbol === g.symbol || state.openGenes.has(g.symbol);
      return `
        <div class="sidebar-gene-item ${isSel ? "is-selected" : ""}" data-sidebar-gene="${g.symbol}" title="${g.symbol} — ${g.name || ""} (${g.variants.length} variants)">
          <span class="sidebar-gene-sym">
            <span class="sidebar-tier-dot" style="background:${dotColor};"></span>
            <span>${g.symbol}</span>
          </span>
          <span class="sidebar-gene-meta">${g.variants.length}v</span>
        </div>
      `;
    }).join("");

    listMount.querySelectorAll("[data-sidebar-gene]").forEach((item) => {
      item.addEventListener("click", () => {
        const sym = item.getAttribute("data-sidebar-gene");
        state.selectedGeneSymbol = sym;
        state.openGenes.add(sym);
        renderActiveReport();
        setTimeout(() => {
          const cardEl = $(`gene-card-${sym}`);
          if (cardEl) cardEl.scrollIntoView({ behavior: "smooth", block: "start" });
        }, 40);
      });
    });
  }

  // ---------------------------------------------------------------------------
  // 5. REPORT SELECTOR, TOOLBAR & VIEW DISPATCHER
  // ---------------------------------------------------------------------------
  function initReportSelector() {
    const sel = $("report-selector");
    const topSel = $("topbar-report-pulldown");
    const scopeSel = $("gene-scope-selector");
    const embeddedKeys = ["dashboard_topn", "ai_clinical", "ai_clinical_old", "ai_deep", "ai_deep_old"];

    if (sel) {
      sel.addEventListener("change", (e) => {
        state.currentReport = e.target.value;
        state.categoryFilter = null;
        if (topSel) {
          topSel.value = embeddedKeys.includes(state.currentReport) ? state.currentReport : "";
        }
        renderActiveReport();
      });
    }
    if (topSel) {
      topSel.addEventListener("change", (e) => {
        const val = e.target.value;
        if (!val) return;
        state.currentReport = val;
        state.categoryFilter = null;
        if (sel) sel.value = val;
        renderActiveReport();
        const ws = $("explorer-workspace");
        if (ws) ws.scrollIntoView({ behavior: "smooth", block: "start" });
      });
    }
    if (scopeSel) {
      scopeSel.addEventListener("change", (e) => {
        state.geneScope = e.target.value || "all";
        if (state.geneScope !== "all") {
          autoExpandScopedGenes(state.geneScope);
        }
        renderActiveReport();
      });
    }
  }

  function autoExpandScopedGenes(scopeMode) {
    const pData = window.COHORT_DATA[state.currentPatient];
    if (!pData) return;
    state.openGenes.clear();
    state.openVariants.clear();
    const scopedGenes = filterAndRankGenes(pData.genes || [], state.currentReport, scopeMode);
    scopedGenes.forEach((g) => {
      state.openGenes.add(g.symbol);
      g.variants.slice(0, 3).forEach((v, idx) => {
        const isT1 = v.tier === "Tier1" || (v.clinvar || "").toLowerCase().includes("pathogenic");
        if (idx === 0 || isT1) {
          state.openVariants.add(v.id);
        }
      });
    });
  }

  function initSearchFilter() {
    const input = $("search-input");
    if (!input) return;
    input.addEventListener("input", (e) => {
      state.searchTerm = e.target.value.trim().toLowerCase();
      renderActiveReport();
    });
  }

  function initToolbarActions() {
    const expandAllBtn = $("btn-expand-all");
    const expandTier1VarsBtn = $("btn-expand-tier1-vars");
    const collapseAllBtn = $("btn-collapse-all");

    if (expandAllBtn) {
      expandAllBtn.addEventListener("click", () => {
        const pData = window.COHORT_DATA[state.currentPatient];
        if (pData) {
          const processed = filterAndRankGenes(pData.genes || [], state.currentReport);
          processed.forEach((g) => state.openGenes.add(g.symbol));
          renderActiveReport();
        }
      });
    }

    if (expandTier1VarsBtn) {
      expandTier1VarsBtn.addEventListener("click", () => {
        const pData = window.COHORT_DATA[state.currentPatient];
        if (!pData) return;
        const processed = filterAndRankGenes(pData.genes || [], state.currentReport);
        let openedCount = 0;
        processed.forEach((g) => {
          const tier1Vars = g.variants.filter((v) =>
            v.tier === "Tier1" ||
            (v.clinvar || "").toLowerCase().includes("pathogenic") ||
            (Number(v.aviPhred) || 0) >= 20 ||
            (Number(v.cadd) || 0) >= 25
          );
          if (tier1Vars.length > 0) {
            state.openGenes.add(g.symbol);
            tier1Vars.forEach((v) => {
              state.openVariants.add(v.id);
              openedCount++;
            });
          }
        });
        if (openedCount === 0) {
          processed.slice(0, 15).forEach((g) => {
            state.openGenes.add(g.symbol);
            g.variants.slice(0, 3).forEach((v) => state.openVariants.add(v.id));
          });
        }
        renderActiveReport();
      });
    }

    if (collapseAllBtn) {
      collapseAllBtn.addEventListener("click", () => {
        state.openGenes.clear();
        state.openVariants.clear();
        state.selectedGeneSymbol = null;
        renderActiveReport();
      });
    }
  }

  function renderActiveReport() {
    const reportType = state.currentReport;
    const container = $("explorer-content-mount");
    const ws = $("explorer-workspace");
    if (!container) return;

    if (["dashboard_topn", "ai_deep", "ai_deep_old", "ai_clinical", "ai_clinical_old"].includes(reportType)) {
      if (ws) ws.classList.add("sidebar-collapsed");
      renderAiSummaryReport(reportType, container);
      return;
    }

    if (ws && !state.sidebarCollapsed) {
      ws.classList.remove("sidebar-collapsed");
    }
    renderScoreReport(reportType, container);
  }

  // ---------------------------------------------------------------------------
  // 6. EMBEDDED AI SUMMARY & DASHBOARD TOP-N REPORTS
  // ---------------------------------------------------------------------------
  function resolveReportRelativeUrl(rawRelPath) {
    if (!rawRelPath) return "";
    const loc = decodeURIComponent(window.location.pathname || "");
    const parts = rawRelPath.split("/");
    const baseName = parts[parts.length - 1];
    if (
      /\/reports\/[^/]+\/[^/]+\.html$/i.test(loc) ||
      loc.includes("Genomics Reports") ||
      loc.includes("/Ontology/")
    ) {
      return baseName;
    }
    if (/\/reports\/[^/]+\.html$/i.test(loc) && rawRelPath.startsWith("reports/")) {
      return rawRelPath.slice("reports/".length);
    }
    return rawRelPath;
  }

  function renderAiSummaryReport(reportType, container) {
    const pData = window.COHORT_DATA[state.currentPatient];
    const isDeep = reportType === "ai_deep" || reportType === "ai_deep_old";
    const isOld = reportType === "ai_deep_old" || reportType === "ai_clinical_old";
    const isTopN = reportType === "dashboard_topn";
    const reportInfo = pData.reportStats[reportType] || { title: "Genomic Report" };

    let rawHtml = pData.meta.clinicalHtml;
    let rawPdf = pData.meta.clinicalPdf || (pData.meta.clinicalHtml || "").replace(/\.html$/i, ".pdf");
    let rawMhtml = "";
    if (isTopN) {
      rawHtml = pData.meta.dashboardTopnHtml || (pData.meta.clinicalHtml || "").replace("_clinical_brief.html", "_dashboard_tier1_topn_report.html");
      rawPdf = pData.meta.dashboardTopnPdf || rawHtml.replace(/\.html$/i, ".pdf");
      rawMhtml = pData.meta.dashboardTopnMhtml || rawHtml.replace(/\.html$/i, ".mhtml");
    } else if (reportType === "ai_clinical_old") {
      rawHtml = pData.meta.clinicalOldHtml || (pData.meta.clinicalHtml || "").replace("_clinical_brief.html", "_clinical_brief (old).html");
      rawPdf = pData.meta.clinicalOldPdf || rawHtml.replace(/\.html$/i, ".pdf");
    } else if (reportType === "ai_deep") {
      rawHtml = pData.meta.deepHtml;
      rawPdf = pData.meta.deepPdf || (pData.meta.deepHtml || "").replace(/\.html$/i, ".pdf");
    } else if (reportType === "ai_deep_old") {
      rawHtml = pData.meta.deepOldHtml || (pData.meta.deepHtml || "").replace("_deep_research_report.html", "_deep_research_report (old).html");
      rawPdf = pData.meta.deepOldPdf || rawHtml.replace(/\.html$/i, ".pdf");
    }

    const reportSrc = resolveReportRelativeUrl(rawHtml);
    const pdfSrc = resolveReportRelativeUrl(rawPdf);
    const mhtmlSrc = rawMhtml ? resolveReportRelativeUrl(rawMhtml) : "";
    const altReportType = isDeep ? (isOld ? "ai_deep" : "ai_deep_old") : (isOld ? "ai_clinical" : "ai_clinical_old");
    const altBtnLabel = isOld ? "Switch to New Edition" : "Switch to (old) Edition";

    let titleLabel = "📋 Executive AI 4-Page Clinic Brief (New Landscape Edition)";
    if (isTopN) {
      titleLabel = "📊 Dashboard Tier-1 &amp; Top-N Multi-View Digest (All 8 Dashboard Report Views)";
    } else if (reportType === "ai_clinical_old") {
      titleLabel = "📋 Executive AI 2-Page Clinic Brief (old — Dual-AI Pathology &amp; Prescribing Guardrails)";
    } else if (reportType === "ai_deep") {
      titleLabel = "📄 Comprehensive AI Deep Research Genomic Dossier (New Edition)";
    } else if (reportType === "ai_deep_old") {
      titleLabel = "📄 AI Deep Research Genomic Dossier (old — Determinative Edition)";
    }

    container.innerHTML = `
      <div class="ai-report-view">
        <div class="ai-report-toolbar">
          <div class="ai-report-toolbar__title">
            ${titleLabel}
            <span style="font-weight:400;color:var(--slate);font-size:12px;margin-left:8px;">Cohort: ${pData.patientName}</span>
          </div>
          <div class="ai-report-toolbar__actions" style="display:flex;gap:6px;flex-wrap:wrap;">
            ${!isTopN ? `<button class="btn-secondary" id="btn-switch-ai-edition" style="border-color:#d97706;color:#b45309;background:#fffbeb;">⇄ ${altBtnLabel}</button>` : ""}
            <a href="${reportSrc}" target="_blank" rel="noopener" class="btn-secondary">Open HTML in New Tab ↗</a>
            ${mhtmlSrc ? `<a href="${mhtmlSrc}" download class="btn-secondary" style="border-color:#7c3aed;color:#6d28d9;background:#f5f3ff;">📦 Download .MHTML</a>` : ""}
            <button class="btn-secondary" onclick="window.frames['ai-frame'].print()">🖨️ Print Report</button>
            <a href="${pdfSrc}" target="_blank" download class="btn-secondary" style="border-color:#059669;color:#047857;background:#ecfdf5;">📥 Save / Open PDF</a>
          </div>
        </div>
        <iframe id="ai-frame" name="ai-frame" class="ai-report-frame" src="${reportSrc}" title="${reportInfo.title}"></iframe>
      </div>
    `;

    const switchBtn = $("btn-switch-ai-edition");
    if (switchBtn) {
      switchBtn.addEventListener("click", () => {
        state.currentReport = altReportType;
        const sel = $("report-selector");
        if (sel) sel.value = altReportType;
        const topSel = $("topbar-report-pulldown");
        if (topSel) topSel.value = altReportType;
        renderActiveReport();
      });
    }
  }

  // ---------------------------------------------------------------------------
  // 7. SCORE-BASED GENE REPORT EXPLORER & TIER-1 / TOP-N SNAPSHOT EXPORTER
  // ---------------------------------------------------------------------------
  function renderScoreReport(reportType, container) {
    const pData = window.COHORT_DATA[state.currentPatient];
    const reportInfo = pData.reportStats[reportType] || { title: "Genomic Explorer", desc: "" };
    const genes = pData.genes || [];

    const allRankedGenes = filterAndRankGenes(genes, reportType, "all");
    const processedGenes = filterAndRankGenes(genes, reportType, state.geneScope);
    renderSidebarGeneNavigator(processedGenes);

    const totalVariantsCount = processedGenes.reduce((acc, g) => acc + g.variants.length, 0);
    const tier1GeneCount = allRankedGenes.filter((g) => g.hasTier1OrPathogenic).length;

    let activeFilterBadges = "";
    if (state.geneScope && state.geneScope !== "all") {
      const scopeLabels = {
        tier1_top15: "🎯 Tier-1 + Top 15 Genes",
        tier1_only: "🔥 Tier-1 / Pathogenic Genes Only",
        top10: "🏆 Top 10 Genes",
        top15: "🏆 Top 15 Genes",
        top25: "🏆 Top 25 Genes"
      };
      activeFilterBadges += `<span class="organ-badge" style="background:#f0fdfa;color:#0f766e;border:1px solid #99f6e4;margin-right:6px;">Scope: <strong>${scopeLabels[state.geneScope] || state.geneScope}</strong> (${processedGenes.length} of ${allRankedGenes.length} genes)</span>`;
    }
    if (state.categoryFilter) {
      activeFilterBadges += `<span class="organ-badge" style="background:#e0f2fe;color:#0369a1;margin-right:6px;">Category: <strong>${state.categoryFilter}</strong></span>`;
    }
    if (state.quadrantFilter) {
      activeFilterBadges += `<span class="organ-badge" style="background:#fff1f2;color:#be123c;margin-right:6px;">Quadrant: <strong>${state.quadrantFilter}</strong></span>`;
    }
    if (state.phaseFilter) {
      activeFilterBadges += `<span class="organ-badge" style="background:#eff6ff;color:#1d4ed8;margin-right:6px;">Phase Audit: <strong>${state.phaseFilter}</strong></span>`;
    }

    const isFocused = state.geneScope !== "all";
    const focusBtnLabel = isFocused
      ? `📚 Show All ${allRankedGenes.length} Genes`
      : `🎯 Focus Tier-1 + Top 15 (${Math.min(Math.max(tier1GeneCount, 15), allRankedGenes.length)} Genes)`;

    let html = `
      <div class="report-header-banner" style="flex-wrap:wrap;gap:12px;">
        <div class="report-header-banner__text" style="flex:1 1 420px;">
          <h2>${reportInfo.title}</h2>
          <p>${reportInfo.desc}</p>
          <div style="margin-top:8px;display:flex;gap:6px;flex-wrap:wrap;align-items:center;">
            <button class="btn-secondary" onclick="window.toggleDashboardTopNScope()" style="font-size:11px;padding:4px 9px;border-color:#0f766e;color:#0f766e;background:#f0fdfa;font-weight:700;">
              ${focusBtnLabel}
            </button>
            <button class="btn-secondary" onclick="window.exportActiveDashboardReportSnapshot('mhtml')" title="Export Tier-1 &amp; Top-N Genes in this Dashboard Report as a self-contained RFC 2557 .MHTML Web Archive" style="font-size:11px;padding:4px 9px;border-color:#7c3aed;color:#6d28d9;background:#f5f3ff;font-weight:700;">
              📦 Export .MHTML (Tier-1 / Top N)
            </button>
            <button class="btn-secondary" onclick="window.exportActiveDashboardReportSnapshot('html')" title="Export Tier-1 &amp; Top-N Genes in this Dashboard Report as a standalone .HTML file" style="font-size:11px;padding:4px 9px;border-color:#0284c7;color:#0369a1;background:#f0f9ff;font-weight:700;">
              📄 Export .HTML (Tier-1 / Top N)
            </button>
            <button class="btn-secondary" onclick="window.exportActiveDashboardReportSnapshot('print')" title="Open clean Print / Save-to-PDF view of just the Tier-1 &amp; Top-N Genes in this Dashboard Report" style="font-size:11px;padding:4px 9px;border-color:#059669;color:#047857;background:#ecfdf5;font-weight:700;">
              🖨️ Print / PDF (Tier-1 / Top N)
            </button>
          </div>
          ${activeFilterBadges ? `<div style="margin-top:6px;">${activeFilterBadges} <button onclick="window.clearAllExplorerFilters()" style="border:none;background:none;color:#dc2626;cursor:pointer;font-size:11px;font-weight:700;margin-left:4px;">✕ Reset Scope &amp; Filters</button></div>` : ""}
        </div>
        <div class="report-header-banner__kpi">
          <div class="kpi-item">
            <span class="kpi-item__num">${processedGenes.length}</span>
            <span class="kpi-item__label">${isFocused ? "Scoped Genes" : "Clinical Genes"}</span>
          </div>
          <div class="kpi-item" style="border-left:1px solid var(--line);padding-left:14px;">
            <span class="kpi-item__num">${totalVariantsCount}</span>
            <span class="kpi-item__label">Annotated Variants</span>
          </div>
        </div>
      </div>
      <div class="gene-catalog-container" id="gene-cards-list">
    `;

    if (processedGenes.length === 0) {
      html += `
        <div style="padding:40px;text-align:center;color:var(--slate);background:#ffffff;border-radius:8px;border:1px solid var(--line);">
          <div style="font-size:24px;margin-bottom:8px;">🔍</div>
          <div style="font-weight:700;">No variants matched the active filter criteria.</div>
          <div style="font-size:12px;margin-top:4px;">Try resetting the gene scope, search query, category filter, or wisdom quadrant filter.</div>
        </div>
      `;
    } else {
      processedGenes.forEach((g) => {
        html += renderGeneCard(g, reportType);
      });
    }

    html += `</div>`;
    container.innerHTML = html;

    attachGeneToggleHandlers();
    attachVariantDrawerHandlers();
    resolveEmptyNcbiSummaries(processedGenes);
  }

  window.toggleDashboardTopNScope = function () {
    const scopeSel = $("gene-scope-selector");
    if (state.geneScope === "all") {
      state.geneScope = "tier1_top15";
      if (scopeSel) scopeSel.value = "tier1_top15";
      autoExpandScopedGenes("tier1_top15");
    } else {
      state.geneScope = "all";
      if (scopeSel) scopeSel.value = "all";
    }
    renderActiveReport();
  };

  window.exportActiveDashboardReportSnapshot = function (mode) {
    const pData = window.COHORT_DATA[state.currentPatient];
    if (!pData) return;
    const reportType = state.currentReport;
    const reportInfo = pData.reportStats[reportType] || { title: "Dashboard Report", desc: "" };
    const effectiveScope = state.geneScope === "all" ? "tier1_top15" : state.geneScope;
    const scopedGenes = filterAndRankGenes(pData.genes || [], reportType, effectiveScope);

    // Temporarily ensure all scoped genes and their top/Tier-1 variants are expanded when rendering snapshot HTML
    const prevGenes = new Set(state.openGenes);
    const prevVars = new Set(state.openVariants);
    scopedGenes.forEach((g) => {
      state.openGenes.add(g.symbol);
      g.variants.forEach((v, idx) => {
        const isT1 = v.tier === "Tier1" || (v.clinvar || "").toLowerCase().includes("pathogenic");
        if (idx < 2 || isT1) {
          state.openVariants.add(v.id);
        }
      });
    });

    const cardsHtml = scopedGenes.map((g) => renderGeneCard(g, reportType)).join("\n");

    // Restore interactive state
    state.openGenes = prevGenes;
    state.openVariants = prevVars;

    let inlineCss = "";
    document.querySelectorAll("style").forEach((s) => { inlineCss += s.innerHTML + "\n"; });

    const snapTitle = `${pData.patientName} — ${reportInfo.title} (Tier-1 & Top-${scopedGenes.length} Genes)`;
    const safeFileStem = `${state.currentPatient}_${reportType}_tier1_top${scopedGenes.length}`;
    const snapHtml = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>${snapTitle}</title>
<style>
${inlineCss}
body { background: #f8fafc; padding: 14px; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }
.btn-gene-action { display: none !important; }
@media print { body { background: #ffffff; padding: 0; } .no-print-snap { display: none !important; } }
</style>
</head>
<body>
  <div class="no-print-snap" style="background:#0f172a;color:#fff;padding:8px 14px;border-radius:6px;margin-bottom:12px;display:flex;justify-content:space-between;align-items:center;">
    <div><strong>${snapTitle}</strong> <span style="font-size:11px;color:#94a3b8;margin-left:8px;">Mean Depth: ${pData.meta?.vcfDepthTelemetry?.mean_depth_called_loci || "45x"}</span></div>
    <button onclick="window.print()" style="background:#2563eb;color:#fff;border:none;padding:5px 12px;border-radius:4px;font-weight:700;cursor:pointer;">🖨️ Print / Save PDF</button>
  </div>
  <div class="report-header-banner" style="margin-bottom:12px;">
    <div class="report-header-banner__text">
      <h2>${reportInfo.title} — Tier-1 &amp; Top-${scopedGenes.length} Focused Report</h2>
      <p>Cohort: <strong>${pData.patientName}</strong> (${state.currentPatient}) · Mean Called Depth: <strong>${pData.meta?.vcfDepthTelemetry?.mean_depth_called_loci || "45x"}</strong> · ${reportInfo.desc}</p>
    </div>
    <div class="report-header-banner__kpi">
      <div class="kpi-item"><span class="kpi-item__num">${scopedGenes.length}</span><span class="kpi-item__label">Focused Genes</span></div>
    </div>
  </div>
  <div class="gene-catalog-container">${cardsHtml}</div>
</body>
</html>`;

    if (mode === "print") {
      const w = window.open("", "_blank");
      if (w) {
        w.document.open();
        w.document.write(snapHtml);
        w.document.close();
        setTimeout(() => { w.focus(); w.print(); }, 350);
      }
      return;
    }

    if (mode === "html") {
      const blob = new Blob([snapHtml], { type: "text/html;charset=utf-8" });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = `${safeFileStem}.html`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      return;
    }

    if (mode === "mhtml") {
      const boundary = "----=_NextPart_000_0000_01DB1CA0.8F4E2B10";
      const nowRfc = new Date().toUTCString();
      const locUrl = `file:///${safeFileStem}.html`;
      const mhtmlParts = [
        "From: <Saved by Genomic Command Center>",
        `Snapshot-Content-Location: ${locUrl}`,
        `Subject: ${snapTitle}`,
        `Date: ${nowRfc}`,
        "MIME-Version: 1.0",
        `Content-Type: multipart/related;\r\n\ttype="text/html";\r\n\tboundary="${boundary}"`,
        "",
        `--${boundary}`,
        "Content-Type: text/html; charset=\"utf-8\"",
        "Content-Transfer-Encoding: 8bit",
        `Content-Location: ${locUrl}`,
        "",
        snapHtml,
        `--${boundary}--`,
        ""
      ];
      const mhtmlBlob = new Blob([mhtmlParts.join("\r\n")], { type: "multipart/related" });
      const a = document.createElement("a");
      a.href = URL.createObjectURL(mhtmlBlob);
      a.download = `${safeFileStem}.mhtml`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }
  };

  window.clearAllExplorerFilters = function () {
    state.geneScope = "all";
    state.categoryFilter = null;
    state.quadrantFilter = null;
    state.phaseFilter = null;
    const scopeSel = $("gene-scope-selector");
    if (scopeSel) scopeSel.value = "all";
    const pData = window.COHORT_DATA[state.currentPatient];
    if (pData) renderWisdomSection(pData);
    renderActiveReport();
  };

  function filterAndRankGenes(genes, reportType, overrideScope) {
    const activeScope = overrideScope !== undefined ? overrideScope : (state.geneScope || "all");
    const list = [];

    genes.forEach((g) => {
      let vars = g.variants || [];

      if (state.categoryFilter) {
        const catLow = state.categoryFilter.toLowerCase();
        const orgLow = (g.organSystem || "").toLowerCase();
        const symLow = (g.symbol || "").toLowerCase();
        const nameLow = (g.name || "").toLowerCase();
        const matchCat = orgLow.includes(catLow) || symLow.includes(catLow) || nameLow.includes(catLow);
        if (!matchCat) return;
      }

      if (state.phaseFilter) {
        const pData = window.COHORT_DATA[state.currentPatient];
        const ctItem = (pData?.wisdom?.cisTransAudit || []).find((x) => x.gene === g.symbol);
        const pStat = ctItem?.phaseState || g.phaseAudit?.phaseStatus || "";
        if (!pStat.includes(state.phaseFilter)) return;
      }

      if (state.quadrantFilter) {
        vars = vars.filter((v) => {
          if (v.concordanceQuadrant === state.quadrantFilter) return true;
          const cv = (v.clinvar || "").toLowerCase();
          const isPlp = cv.includes("pathogenic") && !cv.includes("conflicting") && !cv.includes("uncertain");
          const isProt = cv.includes("protective") || (v.category || "").toLowerCase() === "protective";
          const highAi = (Number(v.cadd) || 0) >= 20 || (Number(v.aviPhred) || 0) >= 20 || (Number(v.revel) || 0) >= 0.5 || (Number(v.alphamissense) || 0) >= 0.56 || (Number(v.spliceai) || 0) >= 0.2;
          if (state.quadrantFilter === "q1_concordant_actionable") return isPlp && highAi;
          if (state.quadrantFilter === "q2_ai_rescue_vus") return !isPlp && !isProt && highAi;
          if (state.quadrantFilter === "q4_protective_offset") return isProt;
          if (state.quadrantFilter === "q3_clinical_mod_low_ai") return !isPlp && !isProt && !highAi;
          return true;
        });
      }

      if (reportType === "rare") {
        vars = vars.filter((v) => v.maf !== null && v.maf !== undefined && Number(v.maf) <= 0.01);
      } else if (reportType === "clinvar") {
        vars = vars.filter((v) => v.clinvar && v.clinvar !== "None" && v.clinvar !== "Not reviewed");
      } else if (reportType === "cadd") {
        vars = vars.filter((v) => v.cadd !== null && v.cadd !== undefined && Number(v.cadd) >= 15.0);
      } else if (reportType === "spliceai") {
        vars = vars.filter((v) => v.spliceai !== null && v.spliceai !== undefined && Number(v.spliceai) >= 0.10);
      } else if (reportType === "atlas") {
        vars = vars.filter((v) => v.aviPhred !== null && v.aviPhred !== undefined);
      } else if (reportType === "alphamissense") {
        vars = vars.filter((v) => v.alphamissense !== null && v.alphamissense !== undefined);
      } else if (reportType === "revel") {
        vars = vars.filter((v) => v.revel !== null && v.revel !== undefined && Number(v.revel) >= 0.20);
      }

      if (state.searchTerm) {
        const q = state.searchTerm;
        const geneMatch = (g.symbol || "").toLowerCase().includes(q) || (g.name || "").toLowerCase().includes(q) || (g.organSystem || "").toLowerCase().includes(q);
        const varMatch = vars.some((v) =>
          (v.id || "").toLowerCase().includes(q) ||
          (v.coordinate || "").toLowerCase().includes(q) ||
          (v.cchange || "").toLowerCase().includes(q) ||
          (v.achange || "").toLowerCase().includes(q) ||
          (v.clinvar || "").toLowerCase().includes(q)
        );
        if (!geneMatch && !varMatch) return;
        if (!geneMatch && varMatch) {
          vars = vars.filter((v) =>
            (v.id || "").toLowerCase().includes(q) ||
            (v.coordinate || "").toLowerCase().includes(q) ||
            (v.cchange || "").toLowerCase().includes(q) ||
            (v.achange || "").toLowerCase().includes(q) ||
            (v.clinvar || "").toLowerCase().includes(q)
          );
        }
      }

      if (vars.length === 0) return;

      const hasTier1OrPathogenic = vars.some((v) => {
        const cv = (v.clinvar || "").toLowerCase();
        return v.tier === "Tier1" || (cv.includes("pathogenic") && !cv.includes("conflicting") && !cv.includes("uncertain"));
      });

      let geneScore = 0;
      let topMetricDisplay = "";

      if (reportType === "cadd") {
        const maxCadd = Math.max(...vars.map((v) => Number(v.cadd) || 0));
        geneScore = maxCadd;
        topMetricDisplay = `Max CADD: ${maxCadd.toFixed(1)}`;
      } else if (reportType === "spliceai") {
        const maxSplice = Math.max(...vars.map((v) => Number(v.spliceai) || 0));
        geneScore = maxSplice;
        topMetricDisplay = `SpliceAI Max: ${maxSplice.toFixed(2)}`;
      } else if (reportType === "atlas") {
        const maxAvi = Math.max(...vars.map((v) => Number(v.aviPhred) || 0));
        geneScore = maxAvi;
        topMetricDisplay = `AlphaGenome: ${maxAvi.toFixed(1)}`;
      } else if (reportType === "alphamissense") {
        const maxAm = Math.max(...vars.map((v) => Number(v.alphamissense) || 0));
        geneScore = maxAm;
        topMetricDisplay = `AlphaMissense: ${maxAm.toFixed(3)}`;
      } else if (reportType === "revel") {
        const maxRev = Math.max(...vars.map((v) => Number(v.revel) || 0));
        geneScore = maxRev;
        topMetricDisplay = `REVEL: ${maxRev.toFixed(3)}`;
      } else if (reportType === "rare") {
        const minMaf = Math.min(...vars.map((v) => (v.maf !== null && v.maf !== undefined ? Number(v.maf) : 1.0)));
        geneScore = 1.0 - minMaf;
        topMetricDisplay = `Min gnomAD AF: ${(minMaf * 100).toFixed(3)}%`;
      } else if (reportType === "clinvar") {
        const cvPriority = { "Pathogenic": 10, "Likely pathogenic": 8, "Protective": 7, "Risk factor": 5, "Uncertain significance": 3 };
        const maxCv = Math.max(...vars.map((v) => cvPriority[v.clinvar] || 1));
        geneScore = maxCv;
        const topCvVar = vars.find((v) => cvPriority[v.clinvar] === maxCv);
        topMetricDisplay = topCvVar ? topCvVar.clinvar : "ClinVar";
      } else {
        const tierPri = { "Tier1": 100, "Tier2": 50, "Tier3": 10 };
        const maxTier = Math.max(...vars.map((v) => tierPri[v.tier] || 1));
        const maxCadd = Math.max(...vars.map((v) => Number(v.cadd) || 0));
        geneScore = maxTier * 10 + maxCadd;
        const topTierVar = vars.find((v) => tierPri[v.tier] === maxTier);
        topMetricDisplay = `${topTierVar ? topTierVar.tier : "Tier3"} · CADD ${maxCadd.toFixed(1)}`;
      }

      vars.sort((a, b) => {
        if (reportType === "cadd") return (Number(b.cadd) || 0) - (Number(a.cadd) || 0);
        if (reportType === "spliceai") return (Number(b.spliceai) || 0) - (Number(a.spliceai) || 0);
        if (reportType === "atlas") return (Number(b.aviPhred) || 0) - (Number(a.aviPhred) || 0);
        if (reportType === "alphamissense") return (Number(b.alphamissense) || 0) - (Number(a.alphamissense) || 0);
        if (reportType === "revel") return (Number(b.revel) || 0) - (Number(a.revel) || 0);
        if (reportType === "rare") return (Number(a.maf) || 0) - (Number(b.maf) || 0);
        return (Number(b.cadd) || 0) - (Number(a.cadd) || 0);
      });

      list.push({
        ...g,
        variants: vars,
        geneScore: geneScore,
        hasTier1OrPathogenic: hasTier1OrPathogenic,
        topMetricDisplay: topMetricDisplay
      });
    });

    list.sort((a, b) => b.geneScore - a.geneScore);

    if (activeScope === "tier1_only") {
      return list.filter((g) => g.hasTier1OrPathogenic);
    }
    if (activeScope === "tier1_top15") {
      const out = [];
      const seen = new Set();
      list.forEach((g) => {
        if (g.hasTier1OrPathogenic) {
          out.push(g);
          seen.add(g.symbol);
        }
      });
      list.forEach((g) => {
        if (!seen.has(g.symbol) && out.length < 15) {
          out.push(g);
          seen.add(g.symbol);
        }
      });
      return out;
    }
    if (activeScope === "top10") return list.slice(0, 10);
    if (activeScope === "top15") return list.slice(0, 15);
    if (activeScope === "top25") return list.slice(0, 25);

    return list;
  }

  /**
   * Renders a Gene Card with Header, Phase Audit Badge, Network Button, Per-Gene Variant Expand, and Table
   */
  function renderGeneCard(gene, reportType) {
    const isExpanded = state.openGenes.has(gene.symbol);
    const varCount = gene.variants.length;
    const summaryText = gene.summary && gene.summary.trim().length > 0 ? gene.summary : (state.cachedNcbiSummaries[gene.symbol] || "");

    const allVarsOpen = gene.variants.every((v) => state.openVariants.has(v.id));
    const pStatus = gene.phaseAudit?.phaseStatus || "";
    let phaseAuditBadge = "";
    if (pStatus.includes("Trans")) {
      phaseAuditBadge = `<span class="badge badge--pathogenic" title="Biallelic Trans Compound Heterozygous Candidates">⚖️ Trans Compound Het</span>`;
    } else if (pStatus.includes("Cis")) {
      phaseAuditBadge = `<span class="badge badge--uncertain" title="Uni-allelic Cis Haplotype (1 WT allele intact)">⚖️ Cis Haplotype</span>`;
    }

    let html = `
      <div class="gene-card ${isExpanded ? "is-expanded" : ""}" id="gene-card-${gene.symbol}">
        <div class="gene-card__header" data-gene-symbol="${gene.symbol}">
          <div class="gene-card__left">
            <span class="gene-toggle-glyph">${isExpanded ? "▼" : "▶"}</span>
            <span class="gene-symbol">${gene.symbol}</span>
            <span class="gene-name">${gene.name || ""}</span>
            <span class="gene-locus">${gene.chrom || ""}:${gene.pos ? Number(gene.pos).toLocaleString() : ""}</span>
            <span class="organ-badge">${gene.organSystem || "Multisystem"}</span>
            ${phaseAuditBadge}
          </div>
          <div class="gene-card__right">
            <button class="btn-gene-action btn-gene-action--network" data-gene-network="${gene.symbol}" title="Open Interactive SVG Gene &amp; Pathway Relationship Network">
              🕸️ Gene Network
            </button>
            <button class="btn-gene-action" data-expand-gene-vars="${gene.symbol}" title="Expand or collapse all variant drawers inside ${gene.symbol}">
              ${allVarsOpen && isExpanded ? "Collapse Variants" : `Expand Variants (${varCount})`}
            </button>
            <span class="gene-metric-badge">${gene.topMetricDisplay}</span>
          </div>
        </div>

        <div class="gene-summary-strip">
          <span class="gene-summary-strip__icon">📖</span>
          <div class="gene-summary-text" id="gene-summary-${gene.symbol}">
            ${summaryText ? summaryText : `<span style="color:var(--slate);font-style:italic;">NCBI RefSeq gene summary indexing in progress…</span>`}
            <span style="display:inline-block;margin-left:8px;font-size:11px;">
              ${gene.omimGene && gene.omimGene !== "100000" ? `<a href="https://omim.org/entry/${gene.omimGene}" target="_blank" rel="noopener" style="margin-right:8px;color:#7c3aed;font-weight:600;">OMIM #${gene.omimGene} ↗</a>` : `<a href="https://omim.org/search?search=${gene.symbol}" target="_blank" rel="noopener" style="margin-right:8px;color:#7c3aed;font-weight:600;">OMIM ↗</a>`}
              ${gene.ncbiGeneId && gene.ncbiGeneId !== "0" ? `<a href="https://www.ncbi.nlm.nih.gov/gene/${gene.ncbiGeneId}" target="_blank" rel="noopener" style="margin-right:8px;color:#0f766e;font-weight:600;">NCBI RefSeq #${gene.ncbiGeneId} ↗</a>` : `<a href="https://www.ncbi.nlm.nih.gov/gene/?term=${gene.symbol}" target="_blank" rel="noopener" style="margin-right:8px;color:#0f766e;font-weight:600;">NCBI Gene ↗</a>`}
              <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${gene.symbol}" target="_blank" rel="noopener" style="color:#b45309;font-weight:600;">GeneCards ↗</a>
            </span>
          </div>
        </div>

        ${isExpanded ? renderGeneVariantTable(gene, reportType) : ""}
      </div>
    `;
    return html;
  }

  function renderGeneVariantTable(gene, reportType) {
    let rowsHtml = "";
    gene.variants.forEach((v) => {
      rowsHtml += renderVariantRowAndDrawer(v, gene, reportType);
    });

    return `
      <div class="gene-variant-table-wrap">
        <table class="variant-table">
          <thead>
            <tr>
              <th style="width:28px;"></th>
              <th>Variant / rsID</th>
              <th>Consequence</th>
              <th>Zygosity / VAF</th>
              <th>Phase</th>
              <th>ClinVar</th>
              <th style="background:#f1f5f9;color:var(--teal-dark);">${getActiveMetricHeader(reportType)}</th>
              <th>CADD</th>
              <th>SpliceAI</th>
              <th>AlphaGenome (AVI &amp; Modality)</th>
              <th>REVEL / AM</th>
              <th>Coordinate</th>
            </tr>
          </thead>
          <tbody>
            ${rowsHtml}
          </tbody>
        </table>
      </div>
    `;
  }

  function getActiveMetricHeader(reportType) {
    if (reportType === "cadd") return "CADD (Primary)";
    if (reportType === "spliceai") return "SpliceAI Δ (Primary)";
    if (reportType === "atlas") return "AVI Phred (Primary)";
    if (reportType === "alphamissense") return "AlphaMissense (Primary)";
    if (reportType === "revel") return "REVEL (Primary)";
    if (reportType === "rare") return "gnomAD AF (Primary)";
    if (reportType === "clinvar") return "ClinVar Review";
    return "Action Tier";
  }

  /**
   * Renders Main Variant Row + Expandable 4-Box Drawer (with AlphaGenome 18-Modality Sparkbars & Track Predictions)
   */
  function renderVariantRowAndDrawer(v, gene, reportType) {
    const isOpen = state.openVariants.has(v.id);
    const isRsId = (v.id || "").startsWith("rs");

    const cv = v.clinvar || "None";
    let cvClass = "badge--benign";
    if (cv.toLowerCase().includes("pathogenic")) cvClass = "badge--pathogenic";
    else if (cv.toLowerCase().includes("protective")) cvClass = "badge--protective";
    else if (cv.toLowerCase().includes("uncertain") || cv.toLowerCase().includes("risk")) cvClass = "badge--uncertain";

    const phase = v.phase || "Unphased";
    const phaseClass = phase.includes("Maternal") ? "badge--phase-maternal" : (phase.includes("Paternal") ? "badge--phase-paternal" : "badge--benign");

    const csqStr = Array.isArray(v.consequence) ? v.consequence.join(", ") : (v.consequence || "MIS");

    let activeMetricVal = "";
    if (reportType === "cadd") {
      activeMetricVal = v.cadd !== null && v.cadd !== undefined ? `<span class="metric-highlight">${Number(v.cadd).toFixed(1)}</span>` : "—";
    } else if (reportType === "spliceai") {
      activeMetricVal = v.spliceai !== null && v.spliceai !== undefined ? `<span class="metric-highlight">${Number(v.spliceai).toFixed(2)}</span>` : "—";
    } else if (reportType === "atlas") {
      activeMetricVal = v.aviPhred !== null && v.aviPhred !== undefined ? `<span class="metric-highlight">Q${Number(v.aviPhred).toFixed(1)}</span>` : "—";
    } else if (reportType === "alphamissense") {
      activeMetricVal = v.alphamissense !== null && v.alphamissense !== undefined ? `<span class="metric-highlight">${Number(v.alphamissense).toFixed(3)}</span>` : "—";
    } else if (reportType === "revel") {
      activeMetricVal = v.revel !== null && v.revel !== undefined ? `<span class="metric-highlight">${Number(v.revel).toFixed(3)}</span>` : "—";
    } else if (reportType === "rare") {
      activeMetricVal = v.maf !== null && v.maf !== undefined ? `<span class="metric-highlight">${(Number(v.maf) * 100).toFixed(3)}%</span>` : "Rare / Private";
    } else if (reportType === "clinvar") {
      activeMetricVal = `<span class="badge ${cvClass}">${v.clinvarRev || "criteria provided"}</span>`;
    } else {
      activeMetricVal = `<span class="badge ${v.tier === "Tier1" ? "badge--pathogenic" : "badge--uncertain"}">${v.tier || "Tier3"}</span>`;
    }

    const caddStr = v.cadd !== null && v.cadd !== undefined ? Number(v.cadd).toFixed(1) : "—";
    const spliceStr = v.spliceai !== null && v.spliceai !== undefined ? Number(v.spliceai).toFixed(2) : "—";
    const aviStr = v.aviPhred !== null && v.aviPhred !== undefined
      ? `<div><span style="font-weight:800;color:#0369a1;">Q${Number(v.aviPhred).toFixed(1)}</span> <span style="font-size:10px;color:var(--slate);">(${v.aviModality || "Multi-Track"})</span></div>`
      : "—";
    const revAmStr = `${v.revel !== null && v.revel !== undefined ? Number(v.revel).toFixed(2) : "—"} / ${v.alphamissense !== null && v.alphamissense !== undefined ? Number(v.alphamissense).toFixed(2) : "—"}`;
    const vafSub = v.vaf !== null && v.vaf !== undefined ? `<div style="font-size:10px;color:var(--slate);">VAF: ${Number(v.vaf).toFixed(2)}${v.dp ? ` · DP:${v.dp}` : ""}</div>` : "";

    const mainRow = `
      <tr class="variant-row ${isOpen ? "is-open" : ""}" data-variant-id="${v.id}">
        <td style="text-align:center;font-weight:700;color:var(--slate);">${isOpen ? "▼" : "▶"}</td>
        <td class="mono" style="font-weight:700;color:var(--teal-dark);">
          ${isRsId ? `<a href="https://www.ncbi.nlm.nih.gov/snp/${v.id}" target="_blank" rel="noopener" style="color:var(--teal-dark);">${v.id}</a>` : v.id}
          ${v.cchange ? `<div style="font-size:11px;font-weight:400;color:var(--slate);">${v.cchange}</div>` : ""}
        </td>
        <td><span class="badge badge--consequence">${csqStr}</span></td>
        <td>${v.zygosity || "Het"} <span class="mono" style="font-size:11px;color:var(--slate);">(${v.genotype || ""})</span>${vafSub}</td>
        <td><span class="badge ${phaseClass}">${phase}</span></td>
        <td><span class="badge ${cvClass}">${cv}</span></td>
        <td style="background:#f8fafc;">${activeMetricVal}</td>
        <td class="mono">${caddStr}</td>
        <td class="mono">${spliceStr}</td>
        <td class="mono">${aviStr}</td>
        <td class="mono" style="font-size:11px;">${revAmStr}</td>
        <td class="mono" style="font-size:11px;color:var(--slate);">${v.coordinate || ""}</td>
      </tr>
    `;

    if (!isOpen) return mainRow;

    const sDetails = v.spliceaiDetails || {};
    const agStr = sDetails.ag !== undefined ? Number(sDetails.ag).toFixed(2) : "0.00";
    const alStr = sDetails.al !== undefined ? Number(sDetails.al).toFixed(2) : "0.00";
    const dgStr = sDetails.dg !== undefined ? Number(sDetails.dg).toFixed(2) : "0.00";
    const dlStr = sDetails.dl !== undefined ? Number(sDetails.dl).toFixed(2) : "0.00";

    const litvarUrl = isRsId
      ? `https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=litvar@${encodeURIComponent(v.id)}%23%23&query=${encodeURIComponent(v.id)}`
      : `https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?query=${encodeURIComponent(v.id)}`;
    const biorxivUrl = `https://www.biorxiv.org/search/${encodeURIComponent((gene.symbol || "") + " " + v.id)}`;
    const europePmcUrl = `https://europepmc.org/search?query=${encodeURIComponent((gene.symbol || "") + " " + v.id)}`;

    const pathologies = gene.associatedPathology || [];
    const topFeatures = v.aviTopFeatures || [];
    const trackPredUrl = v.alphagenomeTrackPredUrl || `https://deepmind.google.com/science/alphagenome/atlas/track-predictions?q=${encodeURIComponent(v.coordinate || "")}&m=variant`;

    const drawerRow = `
      <tr class="variant-drawer-row">
        <td colspan="12">
          <div class="variant-drawer-grid">
            <!-- Box 1: Inheritance & Pathology -->
            <div class="variant-drawer-box">
              <h4>Inheritance &amp; Pathology</h4>
              ${pathologies.length ? pathologies.map((p) => `
                <div style="margin-bottom:6px;">
                  <div style="font-weight:600;font-size:11.5px;">${p.name}</div>
                  <div style="margin-top:2px;display:flex;gap:4px;align-items:center;">
                    ${(p.inheritance || []).map((inh) => `<span class="badge" style="background:#f1f5f9;font-size:10px;">${inh}</span>`).join(" ")}
                    ${p.omim ? `<a href="https://omim.org/entry/${p.omim}" target="_blank" rel="noopener" class="mono" style="font-size:10.5px;color:var(--slate);margin-left:4px;">OMIM #${p.omim} ↗</a>` : ""}
                  </div>
                </div>
              `).join("") : `<div style="font-size:11.5px;color:var(--slate);">Complex / multifactorial quantitative trait modifier.</div>`}
              <div style="margin-top:auto;padding-top:6px;border-top:1px dashed var(--line);font-size:10.5px;color:var(--slate);">
                Quadrant: <strong style="color:var(--ink);">${v.concordanceQuadrant || "Polygenic"}</strong>
              </div>
            </div>

            <!-- Box 2: Transcript & HGVS -->
            <div class="variant-drawer-box">
              <h4>Transcript, Phasing &amp; Depth</h4>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">Transcript ID:</span><span class="drawer-metric-val">${v.transcript || "Canonical"}</span></div>
              ${v.cchange ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">cDNA (HGVS):</span><span class="drawer-metric-val">${v.cchange}</span></div>` : ""}
              ${v.achange ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">Protein (HGVS):</span><span class="drawer-metric-val">${v.achange}</span></div>` : ""}
              ${v.vaf !== undefined && v.vaf !== null ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">Allele Fraction (VAF):</span><span class="drawer-metric-val mono">${Number(v.vaf).toFixed(2)}</span></div>` : ""}
              ${v.dp !== undefined && v.dp !== null ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">Locus Read Depth (DP):</span><span class="drawer-metric-val mono">${v.dp}x</span></div>` : ""}
              ${v.vcfGt ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">VCF Phased GT:</span><span class="drawer-metric-val mono">${v.vcfGt} (${phase})</span></div>` : ""}
            </div>

            <!-- Box 3: SpliceAI & ACMG -->
            <div class="variant-drawer-box">
              <h4>SpliceAI &amp; ACMG Evidence</h4>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">SpliceAI Delta:</span><span class="drawer-metric-val mono">AG:${agStr} AL:${alStr} DG:${dgStr} DL:${dlStr}</span></div>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">Review Status:</span><span class="drawer-metric-val">${v.clinvarRev || "criteria provided"}</span></div>
              ${v.acmgPm5 ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PM5:</span><span class="drawer-metric-val" style="color:var(--concern);">Pathogenic missense at codon</span></div>` : ""}
              ${v.acmgPs1 ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PS1:</span><span class="drawer-metric-val" style="color:var(--concern);">Identical amino acid change known</span></div>` : ""}
              ${v.clinvarId ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ClinVar VCV:</span><span class="drawer-metric-val"><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.clinvarId}/" target="_blank" rel="noopener">VCV#${v.clinvarId} ↗</a></span></div>` : ""}
            </div>

            <!-- Box 4: DeepMind AlphaGenome 18-Modality Attribution & Literature -->
            <div class="variant-drawer-box">
              <h4>AlphaGenome 18-Modality Atlas &amp; Literature</h4>
              ${v.aviPhred !== null && v.aviPhred !== undefined ? `
                <div style="margin-bottom:8px;padding:7px 9px;background:#f0f9ff;border:1px solid #bae6fd;border-radius:5px;font-size:11px;">
                  <div style="display:flex;justify-content:space-between;align-items:center;">
                    <span style="font-weight:800;color:#0369a1;">AlphaGenome AVI Phred:</span>
                    <span style="font-weight:800;font-family:monospace;color:#0284c7;font-size:12.5px;">Q${Number(v.aviPhred).toFixed(1)} <span style="font-size:10px;color:#475569;">(Raw: ${v.aviRaw !== null && v.aviRaw !== undefined ? Number(v.aviRaw).toFixed(3) : "—"})</span></span>
                  </div>
                  <div style="font-size:10.5px;color:#334155;margin-top:2px;display:flex;justify-content:space-between;">
                    <span>Rank: <strong>Top ${v.aviPercentile ? Number(v.aviPercentile).toFixed(2) : "0.00"}%</strong> (Quantile: ${v.aviQuantile ? Number(v.aviQuantile).toFixed(4) : "—"})</span>
                  </div>
                  <div style="font-size:10.5px;color:#0f172a;margin-top:2px;">
                    Primary Driver: <strong style="color:#0369a1;">${v.aviModality || "Multi-Track"}</strong>
                  </div>

                  ${topFeatures.length > 0 ? `
                    <div class="avi-sparkbar-list">
                      ${topFeatures.map((tf) => `
                        <div class="avi-sparkbar-row" title="${tf.key}: ${tf.pct}% relative SHAP attribution (weight ${tf.weight})">
                          <span class="avi-sparkbar-lbl">${tf.label}</span>
                          <div class="avi-sparkbar-track"><div class="avi-sparkbar-fill" style="width:${Math.min(tf.pct, 100)}%;"></div></div>
                          <span class="avi-sparkbar-pct">${tf.pct}%</span>
                        </div>
                      `).join("")}
                    </div>
                  ` : ""}

                  <div style="margin-top:6px;display:flex;gap:6px;flex-wrap:wrap;">
                    <a class="btn-track-pred" href="${trackPredUrl}" target="_blank" rel="noopener" title="Launch DeepMind AlphaGenome Ref vs. Alt RNA-seq, Sashimi Splicing &amp; Chromatin Track Predictions">
                      🧬 Ref vs. Alt Track Predictions ↗
                    </a>
                    <a class="ref-link" style="align-self:center;" href="${v.alphagenomeUrl || `https://deepmind.google.com/science/alphagenome/atlas?q=${v.coordinate}&m=variant`}" target="_blank" rel="noopener">Variant Summary ↗</a>
                  </div>
                </div>
              ` : `<div style="font-size:11px;color:var(--slate);margin-bottom:6px;">AlphaGenome evaluated via SpliceAI / CADD ensemble.</div>`}

              <div style="display:flex;flex-wrap:wrap;gap:6px;margin-bottom:6px;">
                <a class="ref-link" href="${v.ucscUrl || "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38"}" target="_blank" rel="noopener">UCSC Browser ↗</a>
                <a class="ref-link" href="${litvarUrl}" target="_blank" rel="noopener">LitVar2 ↗</a>
                <a class="ref-link" href="${biorxivUrl}" target="_blank" rel="noopener">bioRxiv ↗</a>
                <a class="ref-link" href="${europePmcUrl}" target="_blank" rel="noopener">Europe PMC ↗</a>
              </div>

              ${v.studies && v.studies.length ? v.studies.map((s) => `
                <div style="margin-top:4px;border-top:1px dashed var(--line);padding-top:4px;font-size:10.5px;">
                  <div style="font-weight:700;color:var(--teal-dark);">${s.title}</div>
                  <div style="color:var(--slate);">${s.finding}</div>
                  ${s.pmid ? `<div><a href="https://pubmed.ncbi.nlm.nih.gov/${s.pmid}/" target="_blank" rel="noopener" class="ref-link">PubMed PMID: ${s.pmid} ↗</a></div>` : ""}
                </div>
              `).join("") : `<div style="font-size:10.5px;color:var(--slate);">No specific GWAS study associations flagged.</div>`}
            </div>
          </div>
        </td>
      </tr>
    `;

    return mainRow + drawerRow;
  }

  // ---------------------------------------------------------------------------
  // 8. INTERACTIVE SVG GENE RELATIONSHIP NETWORK MODAL
  // ---------------------------------------------------------------------------
  function initNetworkModal() {
    const modal = $("gene-network-modal");
    const closeBtn = $("btn-close-network-modal");
    if (closeBtn && modal) {
      closeBtn.addEventListener("click", () => {
        modal.style.display = "none";
      });
      modal.addEventListener("click", (e) => {
        if (e.target === modal) modal.style.display = "none";
      });
    }
  }

  function openGeneNetworkModal(symbol) {
    const pData = window.COHORT_DATA[state.currentPatient];
    if (!pData) return;
    const gene = (pData.genes || []).find((g) => g.symbol === symbol);
    if (!gene) return;

    const modal = $("gene-network-modal");
    const titleEl = $("network-modal-title");
    const subEl = $("network-modal-subtitle");
    const svgMount = $("network-svg-mount");
    const infoMount = $("network-info-mount");
    if (!modal || !svgMount || !infoMount) return;

    const net = gene.network || { neighbors: [], stringUrl: `https://string-db.org/cgi/network?identifier=${symbol}&species=9606`, reactomeUrl: `https://reactome.org/content/query?q=${symbol}` };
    const partners = net.neighbors || net.nodes || [];

    titleEl.textContent = `🕸️ ${gene.symbol} (${gene.name || "Clinical Locus"}) — Systems Biology & Co-Mutated Pathway Network`;
    subEl.textContent = `Organ System: ${gene.organSystem || "Multisystem"} · ${gene.variants.length} Patient Variants · ${partners.filter((p) => p.inPatientCallset || p.inPatientCohort).length} Co-Mutated Pathway Partners in Active Cohort`;

    // Build Interactive Radial Force SVG
    const W = 560;
    const H = 440;
    const cx = W / 2;
    const cy = H / 2;
    const R = 152;

    let svg = `<svg viewBox="0 0 ${W} ${H}" style="width:100%;height:auto;max-height:440px;" xmlns="http://www.w3.org/2000/svg">`;
    svg += `
      <defs>
        <radialGradient id="hub-grad" cx="50%" cy="50%" r="50%">
          <stop offset="0%" stop-color="#0f766e" />
          <stop offset="100%" stop-color="#123f41" />
        </radialGradient>
      </defs>
    `;

    // Orbit ring
    svg += `<circle cx="${cx}" cy="${cy}" r="${R}" fill="none" stroke="#cbd5e1" stroke-width="1" stroke-dasharray="4,4" />`;

    // Edges & Partner Nodes
    const count = Math.max(partners.length, 1);
    partners.forEach((p, idx) => {
      const angle = -Math.PI / 2 + (idx / count) * 2 * Math.PI;
      const nx = cx + R * Math.cos(angle);
      const ny = cy + R * Math.sin(angle);
      const isCoMut = Boolean(p.inPatientCallset || p.inPatientCohort);
      const edgeColor = isCoMut ? "#dc2626" : "#64748b";
      const nodeFill = isCoMut ? "#fff1f2" : "#ffffff";
      const nodeStroke = isCoMut ? "#b91c1e" : "#0f766e";

      // Edge line
      svg += `<line x1="${cx}" y1="${cy}" x2="${nx}" y2="${ny}" stroke="${edgeColor}" stroke-width="${isCoMut ? "2.4" : "1.4"}" ${isCoMut ? "" : 'stroke-dasharray="3,2"'} />`;

      // Midpoint edge badge
      const mx = (cx + nx) / 2;
      const my = (cy + ny) / 2;
      svg += `<circle cx="${mx}" cy="${my}" r="4" fill="${edgeColor}" />`;

      // Partner node group
      svg += `<g style="cursor:${isCoMut ? "pointer" : "default"};" data-net-partner="${p.symbol}" data-comut="${isCoMut ? "1" : "0"}">`;
      svg += `<circle cx="${nx}" cy="${ny}" r="28" fill="${nodeFill}" stroke="${nodeStroke}" stroke-width="2.2" filter="drop-shadow(0 2px 4px rgba(0,0,0,0.08))" />`;
      svg += `<text x="${nx}" y="${ny - 2}" text-anchor="middle" font-size="10.5" font-weight="800" font-family="monospace" fill="${isCoMut ? "#991b1b" : "#1e293b"}">${p.symbol}</text>`;
      svg += `<text x="${nx}" y="${ny + 10}" text-anchor="middle" font-size="8" font-weight="600" fill="${isCoMut ? "#dc2626" : "#64748b"}">${isCoMut ? "Co-Mutated" : "Interactome"}</text>`;
      svg += `</g>`;
    });

    // Center Hub Node
    svg += `
      <g>
        <circle cx="${cx}" cy="${cy}" r="38" fill="url(#hub-grad)" stroke="#ffffff" stroke-width="3" filter="drop-shadow(0 4px 10px rgba(15,118,110,0.35))" />
        <text x="${cx}" y="${cy - 3}" text-anchor="middle" font-size="13" font-weight="800" font-family="monospace" fill="#ffffff">${gene.symbol}</text>
        <text x="${cx}" y="${cy + 12}" text-anchor="middle" font-size="9.5" font-weight="600" fill="#ccfbf1">${gene.variants.length} ${gene.variants.length === 1 ? "var" : "vars"}</text>
      </g>
    `;

    // Legend inside SVG
    svg += `
      <g transform="translate(14, ${H - 22})">
        <circle cx="6" cy="-4" r="5" fill="#fff1f2" stroke="#b91c1e" stroke-width="2" />
        <text x="16" y="0" font-size="10" font-weight="700" fill="#991b1b">Co-Mutated in Active Cohort (Click node to jump)</text>
        <circle cx="290" cy="-4" r="5" fill="#ffffff" stroke="#0f766e" stroke-width="2" />
        <text x="300" y="0" font-size="10" font-weight="600" fill="#334155">STRING / Reactome Partner</text>
      </g>
    `;

    svg += `</svg>`;
    svgMount.innerHTML = svg;

    // Right side panel details
    infoMount.innerHTML = `
      <div>
        <div style="font-size:11px;font-weight:800;text-transform:uppercase;color:var(--slate);letter-spacing:0.04em;margin-bottom:6px;">
          External Systems Biology &amp; Pathway Atlases
        </div>
        <div style="display:flex;flex-wrap:wrap;gap:6px;">
          <a href="${net.stringUrl}" target="_blank" rel="noopener" class="btn-secondary" style="font-size:11px;color:#0f766e;border-color:#99f6e4;background:#f0fdfa;">STRING-DB PPI Network ↗</a>
          <a href="${net.reactomeUrl}" target="_blank" rel="noopener" class="btn-secondary" style="font-size:11px;color:#1d4ed8;border-color:#bfdbfe;background:#eff6ff;">Reactome Pathways ↗</a>
          <a href="${net.genemaniaUrl || `https://genemania.org/search/homo-sapiens/${gene.symbol}`}" target="_blank" rel="noopener" class="btn-secondary" style="font-size:11px;">GeneMANIA ↗</a>
          <a href="${net.opentargetsUrl || `https://platform.opentargets.org/search?q=${gene.symbol}&page=1`}" target="_blank" rel="noopener" class="btn-secondary" style="font-size:11px;">OpenTargets ↗</a>
        </div>
      </div>

      <div style="border-top:1px solid var(--line);padding-top:10px;">
        <div style="font-size:11px;font-weight:800;text-transform:uppercase;color:var(--slate);letter-spacing:0.04em;margin-bottom:8px;">
          Connected Functional &amp; Epistatic Partners (${partners.length})
        </div>
        <div style="display:flex;flex-direction:column;gap:6px;">
          ${partners.map((p) => {
            const isCoMut = Boolean(p.inPatientCallset || p.inPatientCohort);
            return `
            <div style="padding:8px 10px;border:1px solid ${isCoMut ? "#fecaca" : "var(--line)"};background:${isCoMut ? "#fff1f2" : "var(--panel-muted)"};border-radius:6px;display:flex;justify-content:space-between;align-items:center;gap:8px;">
              <div>
                <div style="font-family:var(--font-mono);font-weight:800;font-size:12px;color:${isCoMut ? "#991b1b" : "var(--teal-dark)"};">
                  ${p.symbol}
                  ${isCoMut ? `<span class="badge badge--pathogenic" style="margin-left:4px;font-size:9.5px;">Co-Mutated (${p.varCount || 1}v · CADD ${Number(p.maxCadd || 0).toFixed(1)})</span>` : ""}
                </div>
                <div style="font-size:10.5px;color:var(--slate);margin-top:2px;">${p.relation}</div>
              </div>
              ${isCoMut ? `<button class="btn-secondary" data-modal-jump="${p.symbol}" style="font-size:10.5px;padding:3px 8px;">Jump to Gene →</button>` : `<a href="https://www.ncbi.nlm.nih.gov/gene/?term=${p.symbol}" target="_blank" rel="noopener" class="ref-link">NCBI ↗</a>`}
            </div>
          `;
          }).join("")}
        </div>
      </div>
    `;

    modal.style.display = "flex";

    // Click handlers inside modal to jump to co-mutated gene
    svgMount.querySelectorAll("[data-net-partner]").forEach((gEl) => {
      gEl.addEventListener("click", () => {
        if (gEl.getAttribute("data-comut") === "1") {
          const targetSym = gEl.getAttribute("data-net-partner");
          modal.style.display = "none";
          jumpToGeneCard(targetSym);
        }
      });
    });

    infoMount.querySelectorAll("[data-modal-jump]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const targetSym = btn.getAttribute("data-modal-jump");
        modal.style.display = "none";
        jumpToGeneCard(targetSym);
      });
    });
  }

  // ---------------------------------------------------------------------------
  // 9. EVENT HANDLERS & NCBI SUMMARY RESOLVER
  // ---------------------------------------------------------------------------
  function attachGeneToggleHandlers() {
    // Per-gene Expand Variants button
    document.querySelectorAll("[data-expand-gene-vars]").forEach((btn) => {
      btn.addEventListener("click", (e) => {
        e.stopPropagation();
        const sym = btn.getAttribute("data-expand-gene-vars");
        const pData = window.COHORT_DATA[state.currentPatient];
        if (!pData) return;
        const gene = (pData.genes || []).find((g) => g.symbol === sym);
        if (!gene) return;

        state.openGenes.add(sym);
        const allOpen = gene.variants.every((v) => state.openVariants.has(v.id));
        if (allOpen) {
          gene.variants.forEach((v) => state.openVariants.delete(v.id));
        } else {
          gene.variants.forEach((v) => state.openVariants.add(v.id));
        }
        renderActiveReport();
      });
    });

    // Per-gene Network Modal button
    document.querySelectorAll("[data-gene-network]").forEach((btn) => {
      btn.addEventListener("click", (e) => {
        e.stopPropagation();
        const sym = btn.getAttribute("data-gene-network");
        openGeneNetworkModal(sym);
      });
    });

    // Card header accordion click
    document.querySelectorAll(".gene-card__header").forEach((header) => {
      header.addEventListener("click", (e) => {
        if (e.target.closest("button")) return;
        const sym = header.getAttribute("data-gene-symbol");
        if (state.openGenes.has(sym)) {
          state.openGenes.delete(sym);
        } else {
          state.openGenes.add(sym);
        }
        renderActiveReport();
      });
    });
  }

  function attachVariantDrawerHandlers() {
    document.querySelectorAll(".variant-row").forEach((row) => {
      row.addEventListener("click", (e) => {
        if (e.target.tagName.toLowerCase() === "a") return;
        const vid = row.getAttribute("data-variant-id");
        if (state.openVariants.has(vid)) {
          state.openVariants.delete(vid);
        } else {
          state.openVariants.add(vid);
        }
        renderActiveReport();
      });
    });
  }

  async function resolveEmptyNcbiSummaries(genes) {
    const emptyGenes = genes.filter((g) => (!g.summary || g.summary.trim().length === 0 || g.summary.includes("encodes an essential clinical protein")) && !state.cachedNcbiSummaries[g.symbol]);
    if (emptyGenes.length === 0) return;

    for (const g of emptyGenes.slice(0, 20)) {
      try {
        let geneId = g.ncbiGeneId;
        if (!geneId || geneId === "0") {
          const searchRes = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${encodeURIComponent(g.symbol)}[sym]+AND+human[orgn]&retmode=json`);
          const searchData = await searchRes.json();
          const idList = searchData?.esearchresult?.idlist || [];
          if (idList.length > 0) geneId = idList[0];
        }

        if (geneId && geneId !== "0") {
          const sumRes = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=${geneId}&retmode=json`);
          const sumData = await sumRes.json();
          const desc = sumData?.result?.[geneId]?.summary || sumData?.result?.[geneId]?.description || "";
          if (desc) {
            state.cachedNcbiSummaries[g.symbol] = desc;
            const el = $(`gene-summary-${g.symbol}`);
            if (el) {
              const omimLink = (g.omimGene && g.omimGene !== "100000") ? `<a href="https://omim.org/entry/${g.omimGene}" target="_blank" rel="noopener" style="margin-right:8px;color:#7c3aed;font-weight:600;">OMIM #${g.omimGene} ↗</a>` : `<a href="https://omim.org/search?search=${g.symbol}" target="_blank" rel="noopener" style="margin-right:8px;color:#7c3aed;font-weight:600;">OMIM ↗</a>`;
              const ncbiLink = `<a href="https://www.ncbi.nlm.nih.gov/gene/${geneId}" target="_blank" rel="noopener" style="margin-right:8px;color:#0f766e;font-weight:600;">NCBI RefSeq #${geneId} ↗</a>`;
              const gcLink = `<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${g.symbol}" target="_blank" rel="noopener" style="color:#b45309;font-weight:600;">GeneCards ↗</a>`;
              el.innerHTML = `${desc} <span style="display:inline-block;margin-left:8px;font-size:11px;">${omimLink} ${ncbiLink} ${gcLink}</span>`;
            }
          }
        }
      } catch (err) {
        console.warn("NCBI fetch error for " + g.symbol, err);
      }
    }
  }

})();
