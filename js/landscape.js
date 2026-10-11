/**
 * LANDSCAPE CLINICAL DASHBOARD — APPLICATION ENGINE
 * Pure SVG responsive charts, dynamic multi-report table explorer,
 * NCBI gene summary resolution, and dual-cohort patient switching.
 */

(function () {
  "use strict";

  // Application State
  const state = {
    currentPatient: "Daniel_Ehrle",
    currentReport: "visual_explorer",
    searchTerm: "",
    categoryFilter: null,
    barScale: "linear", // "linear" or "log10"
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
        loadPatientData(patientId);
      });
    });
  }

  function loadPatientData(patientId) {
    if (!window.COHORT_DATA || !window.COHORT_DATA[patientId]) {
      console.error("Cohort data not available for:", patientId);
      return;
    }
    const pData = window.COHORT_DATA[patientId];

    // Update Header Meta
    const metaEl = $("topbar-meta");
    if (metaEl) {
      metaEl.textContent = `${pData.patientName} · ${pData.meta.totalVariants} variants · ${pData.meta.totalGenes} genes · Phased WGS`;
    }

    // Render 3 Rows of Charts
    renderAllLandscapeCharts(pData.charts);

    // Render Current Report View
    renderActiveReport();
  }

  // ---------------------------------------------------------------------------
  // 2. SVG CHART ENGINE (PANEL A: STACKED BAR | PANELS B & C: PIE WITH CALLOUTS)
  // ---------------------------------------------------------------------------
  function renderAllLandscapeCharts(charts) {
    if (!charts) return;

    // Groupings by Ontology
    if (charts.row1_ontology) {
      renderStackedBar("chart-row1-bar", charts.row1_ontology.bar, "Diagnostic vs Secondary by HPO Domain", "Variants", true);
      renderPieChart("chart-row1-pie-prior", charts.row1_ontology.piePrior, "Presentation prior to RGS");
      renderPieChart("chart-row1-pie-post", charts.row1_ontology.piePost, "Diagnosis Category after RGS");
    }

    // Groupings by Organ System (only if present in DOM)
    if ($("chart-row2-bar") && charts.row2_organ) {
      renderStackedBar("chart-row2-bar", charts.row2_organ.bar, "Actionability by Anatomical Organ System", "Variants", true);
      renderPieChart("chart-row2-pie-prior", charts.row2_organ.piePrior, "Organ Signals prior to Resolution");
      renderPieChart("chart-row2-pie-post", charts.row2_organ.piePost, "Primary Confirmed Organ Attribution");
    }

    // Groupings by Pharmacogenomics
    if (charts.row3_pharma) {
      renderStackedBar("chart-row3-bar", charts.row3_pharma.bar, "Actionability by Therapeutic Drug Class", "Drug Cards", false);
      renderPieChart("chart-row3-pie-prior", charts.row3_pharma.piePrior, "Metabolizer Phenotypes");
      renderPieChart("chart-row3-pie-post", charts.row3_pharma.piePost, "Clinical Action Categories");
    }
  }

  /**
   * Pure SVG Stacked Bar Chart (Panel A)
   * Matches reference image: Black (Secondary/Non-diagnostic), Red (Actionable/Diagnostic)
   * Displays fraction (e.g. 19/52) above bars, tilted category labels at bottom.
   * Supports Linear & Log10 scaling.
   */
  function renderStackedBar(containerId, data, title, yUnit, allowScaleToggle = false) {
    const mount = $(containerId);
    if (!mount) return;
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

    const n = data.length;
    const colW = chartW / n;
    const barW = Math.max(colW * 0.58, 14);

    let html = "";

    // Optional Log10 Toggle Toolbar
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

    // Horizontal Gridlines & Y-Axis Ticks
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

    // Y-Axis line & label
    svg += `<line x1="${left}" y1="${top}" x2="${left}" y2="${top + chartH}" stroke="#17222b" stroke-width="1.5" />`;
    svg += `<text transform="translate(14, ${top + chartH / 2}) rotate(-90)" text-anchor="middle" font-size="10.5" font-weight="600" fill="#475569">${yUnit} ${isLog ? "(Log₁₀)" : ""}</text>`;

    // X-Axis baseline
    svg += `<line x1="${left}" y1="${top + chartH}" x2="${left + chartW}" y2="${top + chartH}" stroke="#17222b" stroke-width="1.5" />`;

    // Bars
    data.forEach((d, i) => {
      const barX = left + i * colW + (colW - barW) / 2;
      const totalH = (toScale(d.total) / maxScaled) * chartH;
      const actH = (toScale(d.actionable) / maxScaled) * chartH;
      const secH = Math.max(totalH - actH, 0);

      const actY = top + chartH - totalH;
      const secY = top + chartH - secH;

      // Fraction text above bar
      svg += `<text x="${barX + barW / 2}" y="${actY - 6}" text-anchor="middle" font-size="10" font-weight="700" fill="#17222b">${d.fraction}</text>`;

      // Stacked Bar: Red portion (top = Actionable)
      if (actH > 0) {
        svg += `<rect x="${barX}" y="${actY}" width="${barW}" height="${actH}" fill="#B3261E" rx="1" style="cursor:pointer;" data-category="${d.category}">
          <title>${d.category}: ${d.actionable} Actionable / ${d.total} Total (${d.ratePct}%)</title>
        </rect>`;
      }

      // Stacked Bar: Dark portion (bottom = Secondary/Non-diagnostic)
      if (secH > 0) {
        svg += `<rect x="${barX}" y="${secY}" width="${barW}" height="${secH}" fill="#17222B" rx="1" style="cursor:pointer;" data-category="${d.category}">
          <title>${d.category}: ${d.secondary} Secondary / ${d.total} Total</title>
        </rect>`;
      }

      // Category label (Rotated)
      const lblX = barX + barW / 2;
      const lblY = top + chartH + 12;
      const shortCat = d.category.length > 14 ? d.category.slice(0, 13) + "…" : d.category;
      svg += `<text transform="translate(${lblX}, ${lblY}) rotate(-35)" text-anchor="end" font-size="10" font-weight="600" fill="#334155" style="cursor:pointer;" data-category="${d.category}">${shortCat}</text>`;
    });

    // Bottom Legend
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

    // Attach click filters
    mount.querySelectorAll("[data-category]").forEach((el) => {
      el.addEventListener("click", () => {
        const cat = el.getAttribute("data-category");
        setCategoryFilter(cat);
      });
    });

    // Attach scale toggle handlers
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

  /**
   * Pure SVG Pie Chart with Collision-Avoidance Callouts & All-Slice Breakdown Strip
   * Matches reference image panels B and C with white callout cards pointing to slices.
   * Guarantees every single slice displays its exact percentage (.%).
   */
  function renderPieChart(containerId, slices, title) {
    const mount = $(containerId);
    if (!mount) return;
    mount.innerHTML = "";

    const W = 380;
    const H = 265;
    const cx = 190;
    const cy = 132;
    const radius = 72;

    const total = slices.reduce((acc, s) => acc + s.count, 0) || 1;

    let svg = `<svg viewBox="0 0 ${W} ${H}" class="chart-svg" xmlns="http://www.w3.org/2000/svg">`;

    let currentAngle = -Math.PI / 2; // Start from 12 o'clock

    // First pass: Draw slices and gather raw callout positions
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
      const pathD = `M ${cx} ${cy} L ${x1} ${y1} A ${radius} ${radius} 0 ${largeArc} 1 ${x2} ${y2} Z`;

      svg += `<path d="${pathD}" fill="${s.color}" stroke="#ffffff" stroke-width="1.2" style="cursor:pointer;" data-label="${s.label}">
        <title>${s.label}: ${s.pct}% (${s.count}/${total})</title>
      </path>`;

      // All slices with count > 0 receive a callout anchor
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

    // Second pass: Collision Avoidance for Callouts (separate left and right)
    const rightSide = rawCallouts.filter((c) => Math.cos(c.midAngle) >= 0);
    const leftSide = rawCallouts.filter((c) => Math.cos(c.midAngle) < 0);

    rightSide.sort((a, b) => a.outerY - b.outerY);
    leftSide.sort((a, b) => a.outerY - b.outerY);

    function relaxVerticalPositions(list, boxH, gap, minY, maxY) {
      if (list.length === 0) return;
      const minStep = boxH + gap;

      // Downward spacing pass
      for (let i = 1; i < list.length; i++) {
        if (list[i].outerY < list[i - 1].outerY + minStep) {
          list[i].outerY = list[i - 1].outerY + minStep;
        }
      }
      // Upward adjustment if bottom overflows
      if (list[list.length - 1].outerY > maxY) {
        list[list.length - 1].outerY = maxY;
        for (let i = list.length - 2; i >= 0; i--) {
          if (list[i].outerY > list[i + 1].outerY - minStep) {
            list[i].outerY = list[i + 1].outerY - minStep;
          }
        }
      }
      // Top boundary clamp
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

    // Third pass: Draw callout lines & cards with exact .%
    adjustedCallouts.forEach((c) => {
      const isRight = Math.cos(c.midAngle) >= 0;
      const elbowX = isRight ? c.outerX + 14 : c.outerX - 14;
      const boxX = isRight ? elbowX + 2 : elbowX - boxW - 2;
      const boxY = c.outerY - boxH / 2;

      // Pointer polyline
      svg += `<polyline points="${c.innerX},${c.innerY} ${c.outerX},${c.outerY} ${elbowX},${c.outerY}" stroke="#64748b" stroke-width="1" fill="none" />`;

      // Callout card
      svg += `<g style="cursor:pointer;" data-label="${c.slice.label}">`;
      svg += `<rect x="${boxX}" y="${boxY}" width="${boxW}" height="${boxH}" fill="#ffffff" stroke="#94a3b8" stroke-width="0.8" rx="3" filter="drop-shadow(0 1px 2px rgba(0,0,0,0.06))" />`;

      const lblText = c.slice.label.length > 13 ? c.slice.label.slice(0, 12) + "…" : c.slice.label;
      svg += `<text x="${boxX + boxW / 2}" y="${boxY + 10}" text-anchor="middle" font-size="8.5" font-weight="700" fill="#1e293b">${lblText}</text>`;
      svg += `<text x="${boxX + boxW / 2}" y="${boxY + 20}" text-anchor="middle" font-size="8.5" font-weight="600" fill="#0f766e">${c.slice.pct}%, <tspan fill="#64748b">${c.slice.count}</tspan></text>`;
      svg += `</g>`;
    });

    svg += `</svg>`;

    // Fourth pass: Append All-Slices Breakdown Strip below SVG (guarantees every slice percentage is readable)
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

    // Attach click filters
    mount.querySelectorAll("[data-label]").forEach((el) => {
      el.addEventListener("click", () => {
        const lbl = el.getAttribute("data-label");
        setCategoryFilter(lbl);
      });
    });
  }

  function setCategoryFilter(category) {
    if (state.categoryFilter === category) {
      state.categoryFilter = null; // Toggle off
    } else {
      state.categoryFilter = category;
    }
    renderActiveReport();
  }

  // ---------------------------------------------------------------------------
  // 3. REPORT SELECTOR & VIEW DISPATCHER
  // ---------------------------------------------------------------------------
  function initReportSelector() {
    const sel = $("report-selector");
    if (!sel) return;
    sel.addEventListener("change", (e) => {
      state.currentReport = e.target.value;
      state.categoryFilter = null;
      renderActiveReport();
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
    const collapseAllBtn = $("btn-collapse-all");

    if (expandAllBtn) {
      expandAllBtn.addEventListener("click", () => {
        const pData = window.COHORT_DATA[state.currentPatient];
        if (pData) {
          pData.genes.forEach((g) => state.openGenes.add(g.symbol));
          renderActiveReport();
        }
      });
    }

    if (collapseAllBtn) {
      collapseAllBtn.addEventListener("click", () => {
        state.openGenes.clear();
        state.openVariants.clear();
        renderActiveReport();
      });
    }
  }

  function renderActiveReport() {
    const reportType = state.currentReport;
    const container = $("explorer-content-mount");
    if (!container) return;

    // Check if AI Summary Report
    if (reportType === "ai_deep" || reportType === "ai_clinical") {
      renderAiSummaryReport(reportType, container);
      return;
    }

    // Otherwise render Gene-Organized Table Explorer
    renderScoreReport(reportType, container);
  }

  // ---------------------------------------------------------------------------
  // 4. EMBEDDED AI SUMMARY REPORTS (ALT B DOSSIER & ALT C CLINICAL BRIEF)
  // ---------------------------------------------------------------------------
  function renderAiSummaryReport(reportType, container) {
    const pData = window.COHORT_DATA[state.currentPatient];
    const isDeep = reportType === "ai_deep";
    const reportInfo = pData.reportStats[reportType];
    const reportSrc = isDeep ? pData.meta.deepHtml : pData.meta.clinicalHtml;

    container.innerHTML = `
      <div class="ai-report-view">
        <div class="ai-report-toolbar">
          <div class="ai-report-toolbar__title">
            ${isDeep ? "📄 Comprehensive AI Deep Research Genomic Dossier (Alt B)" : "📋 Executive AI Clinical Action Brief (Alt C 3-Page)"}
            <span style="font-weight:400;color:var(--slate);font-size:12px;margin-left:8px;">Patient: ${pData.patientName}</span>
          </div>
          <div class="ai-report-toolbar__actions">
            <a href="${reportSrc}" target="_blank" rel="noopener" class="btn-secondary">Open in New Tab ↗</a>
            <button class="btn-secondary" onclick="window.frames['ai-frame'].print()">Print Report</button>
          </div>
        </div>
        <iframe id="ai-frame" name="ai-frame" class="ai-report-frame" src="${reportSrc}" title="${reportInfo.title}"></iframe>
      </div>
    `;
  }

  // ---------------------------------------------------------------------------
  // 5. SCORE-BASED GENE REPORT EXPLORER
  // ---------------------------------------------------------------------------
  function renderScoreReport(reportType, container) {
    const pData = window.COHORT_DATA[state.currentPatient];
    const reportInfo = pData.reportStats[reportType];
    const genes = pData.genes || [];

    // Filter and score genes
    const processedGenes = filterAndRankGenes(genes, reportType);

    // Build Header Banner
    const totalVariantsCount = processedGenes.reduce((acc, g) => acc + g.variants.length, 0);

    let html = `
      <div class="report-header-banner">
        <div class="report-header-banner__text">
          <h2>${reportInfo.title}</h2>
          <p>${reportInfo.desc}</p>
          ${state.categoryFilter ? `<div style="margin-top:6px;"><span class="organ-badge" style="background:#e0f2fe;color:#0369a1;">Filtered by: <strong>${state.categoryFilter}</strong></span> <button onclick="window.clearCategoryFilter()" style="border:none;background:none;color:#dc2626;cursor:pointer;font-size:11px;font-weight:600;margin-left:6px;">✕ Clear</button></div>` : ""}
        </div>
        <div class="report-header-banner__kpi">
          <div class="kpi-item">
            <span class="kpi-item__num">${processedGenes.length}</span>
            <span class="kpi-item__label">Clinical Genes</span>
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
          <div style="font-size:12px;margin-top:4px;">Try clearing the search query or category filter.</div>
        </div>
      `;
    } else {
      processedGenes.forEach((g) => {
        html += renderGeneCard(g, reportType);
      });
    }

    html += `</div>`;
    container.innerHTML = html;

    // Attach card expansion toggles
    attachGeneToggleHandlers();
    attachVariantDrawerHandlers();

    // Trigger asynchronous NCBI description resolver for empty summaries
    resolveEmptyNcbiSummaries(processedGenes);
  }

  window.clearCategoryFilter = function () {
    state.categoryFilter = null;
    renderActiveReport();
  };

  /**
   * Filters and sorts genes by the active report metric
   */
  function filterAndRankGenes(genes, reportType) {
    const list = [];

    genes.forEach((g) => {
      let vars = g.variants || [];

      // Category filter from chart click
      if (state.categoryFilter) {
        const catLow = state.categoryFilter.toLowerCase();
        const orgLow = (g.organSystem || "").toLowerCase();
        const symLow = (g.symbol || "").toLowerCase();
        const nameLow = (g.name || "").toLowerCase();
        const matchCat = orgLow.includes(catLow) || symLow.includes(catLow) || nameLow.includes(catLow);
        if (!matchCat) return;
      }

      // Metric filter based on report
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

      // Search term filter
      if (state.searchTerm) {
        const q = state.searchTerm;
        const geneMatch = (g.symbol || "").toLowerCase().includes(q) || (g.name || "").toLowerCase().includes(q);
        const varMatch = vars.some((v) =>
          (v.id || "").toLowerCase().includes(q) ||
          (v.coordinate || "").toLowerCase().includes(q) ||
          (v.cchange || "").toLowerCase().includes(q) ||
          (v.achange || "").toLowerCase().includes(q)
        );
        if (!geneMatch && !varMatch) return;
        if (!geneMatch && varMatch) {
          vars = vars.filter((v) =>
            (v.id || "").toLowerCase().includes(q) ||
            (v.coordinate || "").toLowerCase().includes(q) ||
            (v.cchange || "").toLowerCase().includes(q) ||
            (v.achange || "").toLowerCase().includes(q)
          );
        }
      }

      if (vars.length === 0) return;

      // Calculate gene impact score for sorting
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
        topMetricDisplay = `AlphaGenome: Q${maxAvi.toFixed(1)}`;
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
        geneScore = 1.0 - minMaf; // Rarest first
        topMetricDisplay = `Min gnomAD AF: ${(minMaf * 100).toFixed(3)}%`;
      } else if (reportType === "clinvar") {
        const cvPriority = { "Pathogenic": 10, "Likely pathogenic": 8, "Protective": 7, "Risk factor": 5, "Uncertain significance": 3 };
        const maxCv = Math.max(...vars.map((v) => cvPriority[v.clinvar] || 1));
        geneScore = maxCv;
        const topCvVar = vars.find((v) => cvPriority[v.clinvar] === maxCv);
        topMetricDisplay = topCvVar ? topCvVar.clinvar : "ClinVar";
      } else {
        // Visual Explorer: Tier priority + CADD
        const tierPri = { "Tier1": 100, "Tier2": 50, "Tier3": 10 };
        const maxTier = Math.max(...vars.map((v) => tierPri[v.tier] || 1));
        const maxCadd = Math.max(...vars.map((v) => Number(v.cadd) || 0));
        geneScore = maxTier * 10 + maxCadd;
        const topTierVar = vars.find((v) => tierPri[v.tier] === maxTier);
        topMetricDisplay = `${topTierVar ? topTierVar.tier : "Tier3"} · CADD ${maxCadd.toFixed(1)}`;
      }

      // Sort variants within gene
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
        topMetricDisplay: topMetricDisplay
      });
    });

    // Sort genes descending by impact score
    list.sort((a, b) => b.geneScore - a.geneScore);
    return list;
  }

  /**
   * Renders a Gene Card with Header, NCBI Summary, and Variant Table
   */
  function renderGeneCard(gene, reportType) {
    const isExpanded = state.openGenes.has(gene.symbol);
    const varCount = gene.variants.length;
    const summaryText = gene.summary && gene.summary.trim().length > 0 ? gene.summary : (state.cachedNcbiSummaries[gene.symbol] || "");

    let html = `
      <div class="gene-card ${isExpanded ? "is-expanded" : ""}" id="gene-card-${gene.symbol}">
        <div class="gene-card__header" data-gene-symbol="${gene.symbol}">
          <div class="gene-card__left">
            <span class="gene-toggle-glyph">${isExpanded ? "▼" : "▶"}</span>
            <span class="gene-symbol">${gene.symbol}</span>
            <span class="gene-name">${gene.name || ""}</span>
            <span class="gene-locus">${gene.chrom || ""}:${gene.pos ? Number(gene.pos).toLocaleString() : ""}</span>
            <span class="organ-badge">${gene.organSystem || "Multisystem"}</span>
          </div>
          <div class="gene-card__right">
            <span class="gene-metric-badge">${gene.topMetricDisplay}</span>
            <span class="gene-var-count">${varCount} ${varCount === 1 ? "variant" : "variants"}</span>
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

  /**
   * Renders the Variant Table under an expanded Gene
   */
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
              <th>Zygosity</th>
              <th>Phase</th>
              <th>ClinVar</th>
              <th style="background:#f1f5f9;color:var(--teal-dark);">${getActiveMetricHeader(reportType)}</th>
              <th>CADD</th>
              <th>SpliceAI</th>
              <th>AlphaGenome</th>
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
   * Renders Main Variant Row + Expandable 4-Box Drawer (Visual Explorer style)
   */
  function renderVariantRowAndDrawer(v, gene, reportType) {
    const isOpen = state.openVariants.has(v.id);
    const isRsId = (v.id || "").startsWith("rs");

    // ClinVar badge class
    const cv = v.clinvar || "None";
    let cvClass = "badge--benign";
    if (cv.toLowerCase().includes("pathogenic")) cvClass = "badge--pathogenic";
    else if (cv.toLowerCase().includes("protective")) cvClass = "badge--protective";
    else if (cv.toLowerCase().includes("uncertain")) cvClass = "badge--uncertain";

    // Phase badge
    const phase = v.phase || "Unphased";
    const phaseClass = phase === "Maternal" ? "badge--phase-maternal" : (phase === "Paternal" ? "badge--phase-paternal" : "badge--benign");

    // Consequence string
    const csqStr = Array.isArray(v.consequence) ? v.consequence.join(", ") : (v.consequence || "MIS");

    // Active Metric Column Value
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

    // Other metric pills
    const caddStr = v.cadd !== null && v.cadd !== undefined ? Number(v.cadd).toFixed(1) : "—";
    const spliceStr = v.spliceai !== null && v.spliceai !== undefined ? Number(v.spliceai).toFixed(2) : "—";
    const aviStr = v.aviPhred !== null && v.aviPhred !== undefined ? `Q${Number(v.aviPhred).toFixed(1)}` : "—";
    const revAmStr = `${v.revel !== null && v.revel !== undefined ? Number(v.revel).toFixed(2) : "—"} / ${v.alphamissense !== null && v.alphamissense !== undefined ? Number(v.alphamissense).toFixed(2) : "—"}`;

    const mainRow = `
      <tr class="variant-row ${isOpen ? "is-open" : ""}" data-variant-id="${v.id}">
        <td style="text-align:center;font-weight:700;color:var(--slate);">${isOpen ? "▼" : "▶"}</td>
        <td class="mono" style="font-weight:700;color:var(--teal-dark);">
          ${isRsId ? `<a href="https://www.ncbi.nlm.nih.gov/snp/${v.id}" target="_blank" rel="noopener" style="color:var(--teal-dark);">${v.id}</a>` : v.id}
          ${v.cchange ? `<div style="font-size:11px;font-weight:400;color:var(--slate);">${v.cchange}</div>` : ""}
        </td>
        <td><span class="badge badge--consequence">${csqStr}</span></td>
        <td>${v.zygosity || "Het"} <span class="mono" style="font-size:11px;color:var(--slate);">(${v.genotype || ""})</span></td>
        <td><span class="badge ${phaseClass}">${phase}</span></td>
        <td><span class="badge ${cvClass}">${cv}</span></td>
        <td style="background:#f8fafc;">${activeMetricVal}</td>
        <td class="mono">${caddStr}</td>
        <td class="mono">${spliceStr}</td>
        <td class="mono" style="font-weight:600;color:#0369a1;">${aviStr}</td>
        <td class="mono" style="font-size:11px;">${revAmStr}</td>
        <td class="mono" style="font-size:11px;color:var(--slate);">${v.coordinate || ""}</td>
      </tr>
    `;

    if (!isOpen) return mainRow;

    // 4-Box Drawer rendering
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

    const drawerRow = `
      <tr class="variant-drawer-row">
        <td colspan="12">
          <div class="variant-drawer-grid">
            <!-- Box 1: Inheritance & Pathology -->
            <div class="variant-drawer-box">
              <h4>Inheritance & Pathology</h4>
              ${pathologies.length ? pathologies.map((p) => `
                <div style="margin-bottom:6px;">
                  <div style="font-weight:600;font-size:11.5px;">${p.name}</div>
                  <div style="margin-top:2px;display:flex;gap:4px;align-items:center;">
                    ${(p.inheritance || []).map((inh) => `<span class="badge" style="background:#f1f5f9;font-size:10px;">${inh}</span>`).join(" ")}
                    ${p.omim ? `<a href="https://omim.org/entry/${p.omim}" target="_blank" rel="noopener" class="mono" style="font-size:10.5px;color:var(--slate);margin-left:4px;">OMIM #${p.omim} ↗</a>` : ""}
                  </div>
                </div>
              `).join("") : `<div style="font-size:11.5px;color:var(--slate);">Complex / multifactorial quantitative trait modifier.</div>`}
            </div>

            <!-- Box 2: Transcript & HGVS -->
            <div class="variant-drawer-box">
              <h4>Transcript & Nomenclature</h4>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">Transcript ID:</span><span class="drawer-metric-val">${v.transcript || "Canonical"}</span></div>
              ${v.cchange ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">cDNA (HGVS):</span><span class="drawer-metric-val">${v.cchange}</span></div>` : ""}
              ${v.achange ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">Protein (HGVS):</span><span class="drawer-metric-val">${v.achange}</span></div>` : ""}
              ${v.vaf !== undefined && v.vaf !== null ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">Allele Fraction (VAF):</span><span class="drawer-metric-val">${Number(v.vaf).toFixed(2)}</span></div>` : ""}
              ${v.vcfGt ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">VCF Phased GT:</span><span class="drawer-metric-val mono">${v.vcfGt} (${phase})</span></div>` : ""}
            </div>

            <!-- Box 3: SpliceAI & ACMG -->
            <div class="variant-drawer-box">
              <h4>SpliceAI & ACMG Evidence</h4>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">SpliceAI Delta:</span><span class="drawer-metric-val">AG:${agStr} AL:${alStr} DG:${dgStr} DL:${dlStr}</span></div>
              <div class="drawer-metric-row"><span class="drawer-metric-lbl">Review Status:</span><span class="drawer-metric-val">${v.clinvarRev || "criteria provided"}</span></div>
              ${v.acmgPm5 ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PM5:</span><span class="drawer-metric-val" style="color:var(--concern);">Pathogenic missense at codon</span></div>` : ""}
              ${v.acmgPs1 ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ACMG PS1:</span><span class="drawer-metric-val" style="color:var(--concern);">Identical amino acid change known</span></div>` : ""}
              ${v.clinvarId ? `<div class="drawer-metric-row"><span class="drawer-metric-lbl">ClinVar VCV:</span><span class="drawer-metric-val"><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.clinvarId}/" target="_blank" rel="noopener">VCV#${v.clinvarId} ↗</a></span></div>` : ""}
            </div>

            <!-- Box 4: AlphaGenome & Literature -->
            <div class="variant-drawer-box">
              <h4>AlphaGenome & Literature</h4>
              ${v.aviPhred !== null && v.aviPhred !== undefined ? `
                <div style="margin-bottom:8px;padding:6px 8px;background:#f0f9ff;border:1px solid #bae6fd;border-radius:4px;font-size:11px;">
                  <div style="display:flex;justify-content:space-between;align-items:center;">
                    <span style="font-weight:700;color:#0369a1;">DeepMind AlphaGenome:</span>
                    <span style="font-weight:800;font-family:monospace;color:#0284c7;">Q${Number(v.aviPhred).toFixed(1)}</span>
                  </div>
                  <div style="font-size:10.5px;color:#334155;margin-top:2px;">
                    Rank: Top ${v.aviPercentile ? Number(v.aviPercentile).toFixed(2) : "0.00"}% · Modality: <span style="font-weight:700;">${v.aviModality || "Functional"}</span>
                  </div>
                </div>
              ` : `<div style="font-size:11px;color:var(--slate);margin-bottom:6px;">AlphaGenome evaluated via SpliceAI / CADD ensemble.</div>`}

              <div style="display:flex;flex-wrap:wrap;gap:6px;margin-bottom:6px;">
                <a class="ref-link" href="${v.ucscUrl || "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38"}" target="_blank" rel="noopener">UCSC Browser ↗</a>
                <a class="ref-link" href="${v.alphagenomeUrl || `https://deepmind.google.com/science/alphagenome/atlas?q=${v.coordinate}&m=variant`}" target="_blank" rel="noopener">AlphaGenome Atlas ↗</a>
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
  // 6. EVENT HANDLERS & NCBI SUMMARY RESOLVER
  // ---------------------------------------------------------------------------
  function attachGeneToggleHandlers() {
    document.querySelectorAll(".gene-card__header").forEach((header) => {
      header.addEventListener("click", () => {
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
        // Prevent toggle if clicking direct link
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

  /**
   * Asynchronously queries NCBI Entrez e-utilities if a gene description is missing
   */
  async function resolveEmptyNcbiSummaries(genes) {
    const emptyGenes = genes.filter((g) => (!g.summary || g.summary.trim().length === 0 || g.summary.includes("encodes an essential clinical protein")) && !state.cachedNcbiSummaries[g.symbol]);
    if (emptyGenes.length === 0) return;

    for (const g of emptyGenes.slice(0, 20)) { // Batch limit
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
