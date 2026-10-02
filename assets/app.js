/* Portale monitoraggio CavaSA — visualizzazione e confronto di serie temporali.
 * Legge catalog.json e i CSV in data/ (formato standard, orari in UTC+1).
 * Nessun backend: tutto avviene nel browser.
 */
(() => {
  "use strict";

  // ---------------------------------------------------------------- costanti
  const MAX_SERIES = 12;
  const DAY = 86400000;
  const VAR_ORDER = ["rain", "rain_cum", "stage", "vwc", "sm", "psi", "t_soil"];
  const PANEL_WEIGHT = { rain: 0.6 };
  const AUTO_RES = (spanDays) => (spanDays > 90 ? "1d" : spanDays > 4 ? "1h" : "raw");
  const RAW_MAX_DAYS = 62;
  const RES_LABEL = { raw: "dato originale", "1h": "oraria", "1d": "giornaliera" };
  const SHORT = { rain: "Pioggia", rain_cum: "Pioggia cumulata", stage: "Livello idrometrico", vwc: "Contenuto d'acqua", sm: "Umidità terreno (Aranet)",
                  psi: "Potenziale matriciale", t_soil: "Temperatura terreno" };
  const UNIT = (u) => u.replace("m3/m3", "m³/m³");
  const DEFAULT = {
    selected: ["CF_18925_RAIN", "SAL_NV_VWC_015", "SAL_NV_VWC_060", "SAL_NV_PSI_030"],
    from: "2025-10-15 00:00", to: "2025-12-01 00:00",
  };

  // ---------------------------------------------------------------- stato
  const state = {
    selected: [], colors: {}, from: null, to: null,
    res: "auto", mode: "panels", logSuction: false, splitRain: false,
    fStation: "", fVariable: "",
  };
  let catalog, seriesById = {}, stationById = {};
  const cache = new Map();          // url -> Promise<parsed csv>
  let loaded = {};                  // sid -> {res, t[], ts[], v[], f[], min[], max[], cov[]}
  let ignoreRelayout = 0, relayoutTimer = null, renderToken = 0;

  const $ = (id) => document.getElementById(id);
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  // ---------------------------------------------------------------- tempo
  // Timestamp naïve "YYYY-MM-DD HH:MM" (UTC+1) <-> millisecondi trattati come UTC:
  // nessuna conversione di fuso, l'ora mostrata è quella del dato.
  function tms(s) {
    return Date.UTC(+s.slice(0, 4), +s.slice(5, 7) - 1, +s.slice(8, 10),
      +(s.slice(11, 13) || 0), +(s.slice(14, 16) || 0), +(s.slice(17, 19) || 0));
  }
  function tstr(ms) {
    return new Date(ms).toISOString().slice(0, 16).replace("T", " ");
  }
  const fmtDate = (ms) => new Date(ms).toISOString().slice(0, 10);

  // ---------------------------------------------------------------- dati
  function parseCSV(text) {
    const lines = text.split("\n");
    const head = lines[0].trim().split(",");
    const cols = Object.fromEntries(head.map((h) => [h, []]));
    for (let i = 1; i < lines.length; i++) {
      const line = lines[i];
      if (!line) continue;
      const parts = line.split(",");
      for (let j = 0; j < head.length; j++) {
        const h = head[j];
        cols[h].push(h === "time" ? parts[j] : parts[j] === "" ? null : +parts[j]);
      }
    }
    cols.t = cols.time.map(tms);
    return cols;
  }
  function fetchCSV(url) {
    if (!cache.has(url)) {
      cache.set(url, fetch(url).then((r) => {
        if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
        return r.text();
      }).then(parseCSV).catch((e) => { cache.delete(url); throw e; }));
    }
    return cache.get(url);
  }

  async function loadSeries(s, res, t0, t1) {
    let parts;
    if (res === "raw") {
      const y0 = new Date(t0).getUTCFullYear(), y1 = new Date(t1).getUTCFullYear();
      const years = s.data.years.filter((y) => y >= y0 && y <= y1);
      parts = await Promise.all(years.map((y) => fetchCSV(s.data.raw.replace("{year}", y))));
    } else {
      parts = [await fetchCSV(res === "1h" ? s.data.agg_1h : s.data.agg_1d)];
    }
    const out = { res, t: [], ts: [], v: [], f: [], min: [], max: [], cov: [] };
    for (const p of parts) {
      for (let i = 0; i < p.t.length; i++) {
        const t = p.t[i];
        if (t < t0 || t > t1) continue;
        if (res === "raw" && p.flag[i] >= 2) continue;     // dati errati esclusi
        out.t.push(t); out.ts.push(p.time[i]); out.v.push(p.value[i]);
        out.f.push(p.flag ? p.flag[i] : 0);      // aggregati: presente solo per le serie derivate (cumulata)
        out.min.push(p.min ? p.min[i] : null);
        out.max.push(p.max ? p.max[i] : null);
        out.cov.push(p.coverage ? p.coverage[i] : null);
      }
    }
    return out;
  }

  function stepMs(s, res) {
    return res === "raw" ? s.step_min * 60000 : res === "1h" ? 3600000 : DAY;
  }

  // ---------------------------------------------------------------- risoluzione
  function currentRes() {
    const span = (tms(state.to) - tms(state.from)) / DAY;
    if (state.res === "auto") return { res: AUTO_RES(span), note: "" };
    if (state.res === "raw" && span > RAW_MAX_DAYS)
      return { res: "1h", note: `Dato originale disponibile per intervalli fino a ${RAW_MAX_DAYS} giorni: uso la risoluzione oraria.` };
    return { res: state.res, note: "" };
  }

  // ---------------------------------------------------------------- colori
  function assignColor(sid) {
    if (state.colors[sid]) return;
    const used = new Set(state.selected.map((x) => state.colors[x]));
    for (let k = 1; k <= MAX_SERIES; k++) if (!used.has(k)) { state.colors[sid] = k; return; }
  }
  const colorOf = (sid) => css(`--s${state.colors[sid] || 1}`);

  // ---------------------------------------------------------------- etichette
  function axisTitle(variable, short) {
    const v = catalog.variables[variable];
    if (variable === "psi" && state.logSuction) return short ? "kPa (log)" : "Suzione |ψ| (kPa, log)";
    return short ? UNIT(v.unit) : `${SHORT[variable]} (${UNIT(v.unit)})`;
  }
  const panelTitle = (variable) => variable === "psi" && state.logSuction ? "Suzione |ψ|" : catalog.variables[variable].label_it;
  function shortLabel(s) {
    const st = stationById[s.station_id];
    const d = s.depth_m != null ? ` ${Math.round(s.depth_m * 100)} cm` : "";
    return `${SHORT[s.variable]}${d} · ${st.name.replace(/^Salerno\s*–\s*/, "")}`;
  }
  const fmtNum = (x, variable) => x == null || !isFinite(x) ? "–"
    : variable === "vwc" ? x.toFixed(3) : Math.abs(x) >= 100 ? x.toFixed(0) : x.toFixed(1);

  // ---------------------------------------------------------------- gruppi (pannelli / assi)
  // Un gruppo = un asse Y. Di norma uno per grandezza; in modalità Pannelli, con
  // "Piogge in pannelli separati", ogni serie di pioggia ha il proprio pannello.
  function makeGroups() {
    const out = [];
    for (const v of VAR_ORDER) {
      const sids = state.selected.filter((x) => seriesById[x].variable === v);
      if (!sids.length) continue;
      if (v === "rain" && state.splitRain && state.mode === "panels" && sids.length > 1) {
        sids.forEach((sid) => out.push({
          key: `rain:${sid}`, variable: v, sids: [sid],
          title: `Pioggia · ${stationById[seriesById[sid].station_id].name.replace(/^Salerno\s*–\s*/, "")}`,
        }));
      } else {
        out.push({ key: v, variable: v, sids, title: panelTitle(v) });
      }
    }
    return out;
  }

  // ---------------------------------------------------------------- tracce
  function buildTraces(groups, axisOf) {
    const traces = [];
    for (const g of groups) {
      const variable = g.variable;
      for (const sid of g.sids) {
        const s = seriesById[sid], d = loaded[sid];
        if (!d) continue;
        const color = colorOf(sid);
        const ax = axisOf[g.key];
        const transform = variable === "psi" && state.logSuction
          ? (v) => (v != null && v < 0 ? -v : null) : (v) => v;
        const unit = s.unit.replace("m3/m3", "m³/m³");
        const custom = d.t.map((_, i) => [
          fmtNum(d.min[i], variable), fmtNum(d.max[i], variable),
          d.cov[i] == null ? "" : `${Math.round(d.cov[i] * 100)}%`,
          d.f[i] === 1 ? " · cumulata incompleta (dati mancanti)" : ""]);
        const aggInfo = d.res === "raw" ? "" :
          variable === "rain" ? " · max %{customdata[1]} · copertura %{customdata[2]}"
            : variable === "rain_cum" ? " · copertura %{customdata[2]}"
              : " · min %{customdata[0]} max %{customdata[1]} · copertura %{customdata[2]}";
        const name = shortLabel(s);

        if (variable === "rain") {
          traces.push({
            type: "bar", name, x: d.ts, y: d.v, customdata: custom,
            xaxis: "x", yaxis: ax, width: stepMs(s, d.res) * 0.85,
            marker: { color, line: { width: 0 } },
            offset: d.res === "raw" ? -stepMs(s, d.res) * 0.85 : 0,
            hovertemplate: `%{y:.1f} ${unit}${aggInfo}<extra>${name}</extra>`,
          });
        } else {
          // interruzione della linea nelle lacune (> 3 passi)
          const x = [], y = [], cd = [], inc = [], gap = 3 * stepMs(s, d.res);
          for (let i = 0; i < d.t.length; i++) {
            if (i > 0 && d.t[i] - d.t[i - 1] > gap) { x.push(tstr(d.t[i - 1] + 1)); y.push(null); cd.push(["", "", "", ""]); inc.push(false); }
            x.push(d.ts[i]); y.push(transform(d.v[i])); cd.push(custom[i]); inc.push(d.f[i] === 1);
          }
          traces.push({
            type: x.length > 4000 ? "scattergl" : "scatter", mode: "lines", name, x, y, customdata: cd,
            xaxis: "x", yaxis: ax, connectgaps: false,
            line: { color, width: 2 },
            hovertemplate: variable === "rain_cum"
              ? `%{y:.1f} ${unit}${aggInfo}%{customdata[3]}<extra>${name}</extra>`
              : `%{y:.3~f} ${unit}${aggInfo}<extra>${name}</extra>`,
          });
          if (variable === "rain_cum") {
            // tratti con cumulata incompleta (flag 1): tratteggio, ottenuto sovrapponendo
            // alla linea piena una linea tratteggiata del colore dello sfondo
            const ox = [], oy = [];
            let last = -2;
            for (let k = 1; k < x.length; k++) {
              if (!inc[k] || y[k] == null || y[k - 1] == null) continue;
              if (last !== k - 1) { if (ox.length) { ox.push(null); oy.push(null); } ox.push(x[k - 1]); oy.push(y[k - 1]); }
              ox.push(x[k]); oy.push(y[k]); last = k;
            }
            if (ox.length) traces.push({
              type: "scatter", mode: "lines", x: ox, y: oy, xaxis: "x", yaxis: ax, connectgaps: false,
              showlegend: false, hoverinfo: "skip",
              line: { color: css("--surface"), width: 3, dash: "6px,6px" },
            });
          }
          if (d.res === "raw" && variable !== "rain_cum") {
            const sx = [], sy = [];
            d.f.forEach((f, i) => { if (f === 1) { sx.push(d.ts[i]); sy.push(transform(d.v[i])); } });
            if (sx.length) traces.push({
              type: "scatter", mode: "markers", name: `${name} – sospetti`, x: sx, y: sy,
              xaxis: "x", yaxis: ax, showlegend: false,
              marker: { size: 9, color: "rgba(0,0,0,0)", line: { color: css("--suspect"), width: 2 } },
              hovertemplate: `%{y} ${unit} · dato sospetto<extra>${name}</extra>`,
            });
          }
        }
      }
    }
    return traces;
  }

  function baseAxis() {
    return {
      gridcolor: css("--grid"), zeroline: false, linecolor: css("--border"),
      tickfont: { color: css("--text-2"), size: 11 },
      title: { font: { color: css("--text-2"), size: 12 } },
      automargin: true,
    };
  }

  function buildLayout(groups) {
    const layout = {
      paper_bgcolor: css("--surface"), plot_bgcolor: css("--surface"),
      font: { family: "system-ui, -apple-system, Segoe UI, Roboto, sans-serif", color: css("--text") },
      margin: { l: 64, r: 24, t: 40, b: 40 },
      showlegend: true,
      legend: { orientation: "h", x: 0, y: 1.0, yanchor: "bottom", font: { size: 12, color: css("--text-2") } },
      hovermode: "x unified",
      hoverlabel: { bgcolor: css("--surface"), bordercolor: css("--border"), font: { color: css("--text"), size: 12 } },
      bargap: 0, barmode: "overlay", annotations: [],
      xaxis: {
        ...baseAxis(), type: "date", range: [state.from, state.to],
        showspikes: true, spikemode: "across", spikesnap: "cursor", spikethickness: 1,
        spikecolor: css("--muted"), spikedash: "solid",
        hoverformat: "%d/%m/%Y %H:%M",
      },
    };
    const axisOf = {};

    if (state.mode === "panels") {
      const gap = 0.035;
      const w = groups.map((g) => PANEL_WEIGHT[g.variable] || 1);
      const total = w.reduce((a, b) => a + b, 0);
      const usable = 1 - gap * (groups.length - 1);
      let top = 1, firstRain = null;
      groups.forEach((grp, i) => {
        const g = grp.variable;
        const h = usable * w[i] / total;
        const id = i === 0 ? "y" : `y${i + 1}`;
        const key = i === 0 ? "yaxis" : `yaxis${i + 1}`;
        axisOf[grp.key] = id;
        layout.annotations.push({
          text: grp.title, xref: "paper", yref: "paper", x: 0, y: top, xanchor: "left", yanchor: "top",
          showarrow: false, font: { size: 12, color: css("--text-2") }, bgcolor: css("--surface"), borderpad: 2,
        });
        layout[key] = {
          ...baseAxis(), domain: [Math.max(0, top - h), top],
          title: { ...baseAxis().title, text: axisTitle(g, true) },
          type: g === "psi" && state.logSuction ? "log" : "linear",
          ...(g === "psi" && state.logSuction ? { dtick: 1 } : {}),
          rangemode: g === "rain" || g === "rain_cum" ? "tozero" : "normal",
        };
        // pannelli di pioggia separati: stessa scala, per un confronto corretto
        if (g === "rain") { if (firstRain) layout[key].matches = firstRain; else firstRain = id; }
        top -= h + gap;
      });
      layout.xaxis.anchor = groups.length > 1 ? `y${groups.length}` : "y";
    } else {
      // asse multipli: pioggia a destra come ietogramma rovesciato, le altre alternate
      const others = groups.filter((g) => g.variable !== "rain");
      const left = others.filter((_, i) => i % 2 === 0);
      const right = others.filter((_, i) => i % 2 === 1);
      const rainG = groups.find((g) => g.variable === "rain");
      if (rainG) right.unshift(rainG);
      const step = 0.075;
      const x0 = step * Math.max(0, left.length - 1), x1 = 1 - step * Math.max(0, right.length - 1);
      layout.xaxis.domain = [x0, x1];
      layout.margin.r = 64;
      let n = 0;
      const place = (grp, side, k) => {
        const g = grp.variable;
        n += 1;
        const id = n === 1 ? "y" : `y${n}`, key = n === 1 ? "yaxis" : `yaxis${n}`;
        axisOf[grp.key] = id;
        const color = colorOf(grp.sids[0]);
        const ax = {
          ...baseAxis(), side, showgrid: n === 1,
          title: { text: axisTitle(g, false), font: { color, size: 12 } },
          tickfont: { color, size: 11 },
          type: g === "psi" && state.logSuction ? "log" : "linear",
          ...(g === "psi" && state.logSuction ? { dtick: 1 } : {}),
        };
        if (n > 1) ax.overlaying = "y";
        if (k > 0) { ax.anchor = "free"; ax.position = side === "left" ? x0 - step * k : x1 + step * k; }
        if (g === "rain") {
          const maxR = Math.max(1, ...state.selected.filter((x) => seriesById[x].variable === "rain")
            .flatMap((x) => (loaded[x] ? loaded[x].v : [])).filter((v) => v != null));
          ax.range = [maxR * 2.6, 0];          // barre appese in alto, occupano ~40% dell'altezza
          ax.showgrid = false;
        }
        layout[key] = ax;
      };
      left.forEach((g, k) => place(g, "left", k));
      right.forEach((g, k) => place(g, "right", k));
    }
    return { layout, axisOf };
  }

  // ---------------------------------------------------------------- rendering
  async function render() {
    const token = ++renderToken;
    const chart = $("chart");
    writeHash();
    syncInputs();
    if (!state.selected.length) {
      Plotly.purge(chart);
      chart.innerHTML = '<div class="empty">Seleziona una o più serie dall\'elenco a sinistra.</div>';
      $("status").textContent = "";
      renderStats();
      return;
    }
    const { res, note } = currentRes();
    const t0 = tms(state.from), t1 = tms(state.to);
    const pad = Math.max((t1 - t0) * 0.5, DAY);
    $("status").textContent = "Caricamento…";
    try {
      const results = await Promise.all(state.selected.map((sid) =>
        loadSeries(seriesById[sid], res, t0 - pad, t1 + pad).then((d) => [sid, d])));
      if (token !== renderToken) return;                 // superato da una richiesta più recente
      loaded = Object.fromEntries(results);
    } catch (e) {
      $("status").textContent = `Errore nel caricamento dei dati: ${e.message}`;
      return;
    }
    const groups = makeGroups();
    // altezza minima del riquadro: in Pannelli ~130 px per pannello; oltre, l'utente lo ridimensiona liberamente
    const minH = state.mode === "panels" ? 110 + groups.length * 130 : 420;
    chart.style.minHeight = `${minH}px`;
    const { layout, axisOf } = buildLayout(groups);
    const traces = buildTraces(groups, axisOf);
    if (chart.querySelector(".empty")) chart.innerHTML = "";
    ignoreRelayout++;
    await Plotly.react(chart, traces, layout, {
      responsive: true, displaylogo: false, scrollZoom: false, locale: "it",
      modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
      toImageButtonOptions: { filename: "monitoraggio_cavasa", scale: 2 },
    });
    setTimeout(() => { ignoreRelayout = Math.max(0, ignoreRelayout - 1); }, 50);
    if (!chart._relayoutBound) {
      chart.on("plotly_relayout", onRelayout);
      chart._relayoutBound = true;
      // ridimensionamento del riquadro (maniglia in basso a destra o cambio di larghezza):
      // i pannelli sono definiti in frazioni dell'altezza, quindi scalano insieme al riquadro
      let rTimer = null;
      new ResizeObserver(() => {
        clearTimeout(rTimer);
        rTimer = setTimeout(() => { if (chart.data) { ignoreRelayout++; Plotly.Plots.resize(chart).finally(() => setTimeout(() => { ignoreRelayout = Math.max(0, ignoreRelayout - 1); }, 50)); } }, 60);
      }).observe(chart);
    }
    const nPts = Object.values(loaded).reduce((a, d) => a + d.t.length, 0);
    $("status").textContent = `Risoluzione ${RES_LABEL[res]} · ${nPts.toLocaleString("it-IT")} valori caricati`
      + (note ? ` · ${note}` : "") + " · Trascina sul grafico per ingrandire, doppio clic per tornare all'intervallo completo.";
    renderStats();
  }

  function onRelayout(ev) {
    if (ignoreRelayout) return;
    let a = ev["xaxis.range[0]"], b = ev["xaxis.range[1]"];
    if (ev["xaxis.range"]) [a, b] = ev["xaxis.range"];
    if (ev["xaxis.autorange"]) { setFullExtent(); }
    else if (a && b) { state.from = tstr(tms(String(a))); state.to = tstr(tms(String(b))); }
    else return;
    clearTimeout(relayoutTimer);
    relayoutTimer = setTimeout(render, 250);
  }

  function setFullExtent() {
    const sel = state.selected.map((x) => seriesById[x]);
    if (!sel.length) return;
    const a = Math.min(...sel.map((s) => tms(s.start.slice(0, 16).replace("T", " "))));
    const b = Math.max(...sel.map((s) => tms(s.end.slice(0, 16).replace("T", " "))));
    state.from = tstr(a); state.to = tstr(b);
  }

  // ---------------------------------------------------------------- statistiche
  function visibleStats(sid) {
    const s = seriesById[sid], d = loaded[sid];
    if (!d) return null;
    const t0 = tms(state.from), t1 = tms(state.to);
    let n = 0, sum = 0, mn = Infinity, mx = -Infinity, covSum = 0, covN = 0;
    for (let i = 0; i < d.t.length; i++) {
      if (d.t[i] < t0 || d.t[i] > t1 || d.v[i] == null) continue;
      n++; sum += d.v[i];
      // pioggia: estremi del cumulato alla risoluzione visualizzata; altre grandezze: estremi dei dati
      const isRain = s.variable === "rain";
      mn = Math.min(mn, isRain ? d.v[i] : (d.min[i] ?? d.v[i]));
      mx = Math.max(mx, isRain ? d.v[i] : (d.max[i] ?? d.v[i]));
      if (d.cov[i] != null) { covSum += d.cov[i]; covN++; }
    }
    const expected = (t1 - t0) / stepMs(s, d.res);
    const cov = d.res === "raw" ? n / expected : covN ? (covSum / expected) : 0;
    return { n, mean: n ? sum / n : null, sum, min: n ? mn : null, max: n ? mx : null,
             cov: Math.min(1, cov), res: d.res };
  }

  function renderStats() {
    const tbl = $("statsTable");
    if (!state.selected.length) { tbl.innerHTML = ""; return; }
    const rows = state.selected.map((sid) => {
      const s = seriesById[sid], st = visibleStats(sid);
      if (!st) return "";
      const v = s.variable;
      return `<tr><td><span class="swatch" style="background:${colorOf(sid)}"></span>${shortLabel(s)}</td>
        <td>${s.unit.replace("m3/m3", "m³/m³")}</td><td>${RES_LABEL[st.res]}</td>
        <td>${st.n.toLocaleString("it-IT")}</td>
        <td>${fmtNum(st.min, v)}</td><td>${v === "rain" || v === "rain_cum" ? "–" : fmtNum(st.mean, v)}</td><td>${fmtNum(st.max, v)}</td>
        <td>${v === "rain" ? st.sum.toFixed(1) : "–"}</td>
        <td>${Math.round(st.cov * 100)}%</td></tr>`;
    }).join("");
    tbl.innerHTML = `<thead><tr><th>Serie</th><th>Unità</th><th>Risoluzione</th><th>N</th>
      <th>Min</th><th>Media</th><th title="Per la pioggia: massimo cumulato alla risoluzione visualizzata (orario, giornaliero o per singola registrazione). Per la pioggia cumulata dal 1 gennaio: massimo raggiunto nell'intervallo">Max</th><th>Totale</th><th>Copertura</th></tr></thead><tbody>${rows}</tbody>`;
  }

  function exportCSV() {
    const t0 = tms(state.from), t1 = tms(state.to);
    const lines = ["series_id,time_utc+1,value,resolution"];
    for (const sid of state.selected) {
      const d = loaded[sid]; if (!d) continue;
      for (let i = 0; i < d.t.length; i++)
        if (d.t[i] >= t0 && d.t[i] <= t1) lines.push(`${sid},${d.ts[i]},${d.v[i] ?? ""},${d.res}`);
    }
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `cavasa_${state.from.slice(0, 10)}_${state.to.slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  // ---------------------------------------------------------------- stazioni
  const isCF = (st) => /^CF_/.test(st.station_id);
  const elevLabel = (st) => st.elevation_m != null ? `${st.elevation_m} m s.l.m.` : "quota n.d.";
  // ordine: stazioni UNISA, poi rete CF per quota crescente
  const stationOrder = (a, b) => (isCF(a) - isCF(b)) || ((a.elevation_m ?? 0) - (b.elevation_m ?? 0));

  // ---------------------------------------------------------------- elenco serie
  function renderSeriesList() {
    const list = $("seriesList");
    const vis = catalog.series.filter((s) =>
      (!state.fStation || s.station_id === state.fStation) &&
      (!state.fVariable || s.variable === state.fVariable));
    const byStation = {};
    vis.forEach((s) => (byStation[s.station_id] ||= []).push(s));
    const ordered = Object.entries(byStation)
      .sort(([a], [b]) => stationOrder(stationById[a], stationById[b]));
    list.innerHTML = ordered.map(([stId, arr]) => `
      <div class="st-group"><div class="st-title">${stationById[stId].name}
        <span class="st-elev">${stationById[stId].elevation_m != null ? `${stationById[stId].elevation_m} m` : ""}</span></div>
      ${arr.sort((a, b) => VAR_ORDER.indexOf(a.variable) - VAR_ORDER.indexOf(b.variable) || (a.depth_m ?? 0) - (b.depth_m ?? 0))
        .map((s) => {
          const on = state.selected.includes(s.series_id);
          const d = s.depth_m != null ? ` ${Math.round(s.depth_m * 100)} cm` : "";
          return `<div class="s-item" title="${s.label}">
            <label><input type="checkbox" data-sid="${s.series_id}" ${on ? "checked" : ""}>
              <span class="swatch" style="background:${on ? colorOf(s.series_id) : "transparent"}"></span>
              <span class="name">${SHORT[s.variable]}${d}</span></label>
            <button class="info-btn" type="button" data-info="${s.series_id}" title="Metadati">ⓘ</button></div>`;
        }).join("")}</div>`).join("") || '<p class="muted small">Nessuna serie con questi filtri.</p>';
  }

  function toggleSeries(sid, on) {
    if (on) {
      if (state.selected.length >= MAX_SERIES) {
        $("status").textContent = `Puoi confrontare al massimo ${MAX_SERIES} serie.`;
        renderSeriesList(); return;
      }
      assignColor(sid);
      state.selected.push(sid);
    } else {
      state.selected = state.selected.filter((x) => x !== sid);
      delete state.colors[sid];
    }
    renderSeriesList();
    render();
  }

  function showInfo(sid) {
    const s = seriesById[sid], st = stationById[s.station_id], q = s.qc_summary || {};
    const li = (arr) => arr.length ? `<ul>${arr.map((x) => `<li>${x}</li>`).join("")}</ul>` : "–";
    $("infoBody").innerHTML = `
      <h3>${s.label}</h3>
      <dl>
        <dt>Codice serie</dt><dd><code>${s.series_id}</code></dd>
        <dt>Stazione</dt><dd>${st.name}${st.external_code ? ` · codice ${st.external_code}` : ""}</dd>
        <dt>Posizione</dt><dd>${st.lat != null ? `${st.lat.toFixed(5)}, ${st.lon.toFixed(5)}` : "n.d."} · ${elevLabel(st)}</dd>
        <dt>Gestore</dt><dd>${st.owner || "–"}</dd>
        <dt>Sensore</dt><dd>${s.sensor_model || "–"}${s.sensor_serial ? ` · matr. ${s.sensor_serial}` : ""}${s.logger_port ? ` · porta ${s.logger_port}` : ""}</dd>
        <dt>Grandezza</dt><dd>${catalog.variables[s.variable].label_it} [${s.unit}]</dd>
        <dt>Passo</dt><dd>${s.step_min} min</dd>
        <dt>Periodo</dt><dd>${s.start.slice(0, 16).replace("T", " ")} → ${s.end.slice(0, 16).replace("T", " ")} (UTC+1)</dd>
        <dt>Valori</dt><dd>${(q.n_total || s.n_values).toLocaleString("it-IT")} · sospetti ${q.n_suspect ?? 0} · errati ${q.n_bad ?? 0}</dd>
        <dt>Lacune &gt; 24 h</dt><dd>${li((s.gaps_over_24h || []).map((g) => `${g.from.slice(0, 16).replace("T", " ")} → ${g.to.slice(0, 16).replace("T", " ")} (${g.hours} h)`))}</dd>
        <dt>Note di qualità</dt><dd>${li(s.qc_notes || [])}</dd>
        ${s.derived_from ? `<dt>Serie derivata da</dt><dd><code>${s.derived_from}</code></dd>` : ""}
        <dt>File sorgente</dt><dd>${li((s.sources || []).map((x) => `${x.file} · ${x.sheet} · ${x.column}`))}</dd>
      </dl>`;
    $("infoDlg").showModal();
  }

  // ---------------------------------------------------------------- mappa
  function initMap() {
    const withXY = catalog.stations.filter((s) => s.lat != null);
    const without = catalog.stations.filter((s) => s.lat == null);
    $("mapNote").textContent = without.length
      ? `Posizione non disponibile: ${without.map((s) => s.name).join(", ")}.` : "";
    if (!window.L || !withXY.length) { $("map").style.display = "none"; return; }
    const map = L.map("map", { scrollWheelZoom: false, attributionControl: true });
    L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
      maxZoom: 18, attribution: "© OpenStreetMap",
    }).addTo(map);
    // stazioni con le stesse coordinate (o molto vicine, es. Meter NV e V): disposte a ventaglio attorno al punto reale,
    // con scarto fisso in pixel (indicativo, solo per renderle selezionabili)
    const groups = {};
    withXY.forEach((s) => (groups[`${s.lat.toFixed(4)},${s.lon.toFixed(4)}`] ||= []).push(s));
    const items = [];
    const pts = withXY.map((s) => {
      const grp = groups[`${s.lat.toFixed(4)},${s.lon.toFixed(4)}`].slice().sort((x, y) => x.station_id.localeCompare(y.station_id));
      const spread = grp.length > 1;
      const m = L.circleMarker([s.lat, s.lon], {
        radius: isCF(s) ? 7 : 8, weight: 2, color: "#ffffff",
        fillColor: isCF(s) ? css("--s2") : css("--s1"), fillOpacity: 1,
      }).addTo(map);
      m.bindTooltip(`<b>${s.name}</b><br>${elevLabel(s)}${spread ? "<br>posizione indicativa (a pochi metri da " + (grp.length - 1) + " altra/e stazione/i)" : ""}<br><span style="opacity:.75">clic per filtrare le serie</span>`);
      m.on("click", () => { state.fStation = s.station_id; $("fStation").value = s.station_id; renderSeriesList(); });
      if (spread) items.push({ m, s, k: grp.indexOf(s), n: grp.length });
      return [s.lat, s.lon];
    });
    const fan = () => items.forEach(({ m, s, k, n }) => {
      const c = map.latLngToLayerPoint([s.lat, s.lon]), a = -Math.PI / 2 + (2 * Math.PI * k) / n;
      m.setLatLng(map.layerPointToLatLng(L.point(c.x + 10 * Math.cos(a), c.y + 10 * Math.sin(a))));
    });
    map.on("zoomend", fan);
    map.fitBounds(pts, { padding: [30, 30], maxZoom: 14 });
    fan();
  }

  // ---------------------------------------------------------------- URL / controlli
  function writeHash() {
    const p = new URLSearchParams({
      s: state.selected.join(","), from: state.from, to: state.to,
      mode: state.mode, res: state.res, log: state.logSuction ? "1" : "0", split: state.splitRain ? "1" : "0",
    });
    history.replaceState(null, "", `#${p.toString()}`);
  }
  function readHash() {
    const p = new URLSearchParams(location.hash.slice(1));
    const sel = (p.get("s") || "").split(",").filter((x) => seriesById[x]);
    state.selected = sel.length || p.has("s") ? sel : DEFAULT.selected.filter((x) => seriesById[x]);
    state.from = p.get("from") || DEFAULT.from;
    state.to = p.get("to") || DEFAULT.to;
    state.mode = p.get("mode") === "overlay" ? "overlay" : "panels";
    state.res = ["auto", "1h", "1d", "raw"].includes(p.get("res")) ? p.get("res") : "auto";
    state.logSuction = p.get("log") === "1";
    state.splitRain = p.get("split") === "1";
    state.selected.forEach(assignColor);
  }
  function syncInputs() {
    $("dFrom").value = state.from.slice(0, 10);
    $("dTo").value = state.to.slice(0, 10);
    $("resSel").value = state.res;
    $("logSuction").checked = state.logSuction;
    $("splitRain").checked = state.splitRain;
    $("splitRain").disabled = state.mode !== "panels";
    $("splitRain").parentElement.classList.toggle("disabled", state.mode !== "panels");
    document.querySelectorAll(".seg button").forEach((b) => b.classList.toggle("on", b.dataset.mode === state.mode));
  }

  function bindControls() {
    const fill = (el, opts) => { el.innerHTML = opts.map(([v, t]) => `<option value="${v}">${t}</option>`).join(""); };
    fill($("fStation"), [["", "Tutte le stazioni"], ...[...catalog.stations].sort(stationOrder).map((s) => [s.station_id, s.name])]);
    fill($("fVariable"), [["", "Tutte le grandezze"],
      ...VAR_ORDER.filter((v) => catalog.series.some((s) => s.variable === v)).map((v) => [v, catalog.variables[v].label_it])]);
    $("fStation").onchange = (e) => { state.fStation = e.target.value; renderSeriesList(); };
    $("fVariable").onchange = (e) => { state.fVariable = e.target.value; renderSeriesList(); };

    $("seriesList").addEventListener("change", (e) => {
      if (e.target.dataset.sid) toggleSeries(e.target.dataset.sid, e.target.checked);
    });
    $("seriesList").addEventListener("click", (e) => {
      const b = e.target.closest("[data-info]");
      if (b) showInfo(b.dataset.info);
    });
    $("clearBtn").onclick = () => { state.selected = []; state.colors = {}; renderSeriesList(); render(); };

    const setRange = () => {
      const a = $("dFrom").value, b = $("dTo").value;
      if (!a || !b || a >= b) return;
      state.from = `${a} 00:00`; state.to = `${b} 00:00`;
      render();
    };
    $("dFrom").onchange = setRange; $("dTo").onchange = setRange;
    document.querySelectorAll(".presets button").forEach((b) => b.onclick = () => {
      if (b.dataset.days === "all") { setFullExtent(); render(); return; }   // tutto il set di dati delle serie selezionate
      // 7 g / 30 g / 1 anno: l'inizio resta la data del box "Da"; cambia solo la fine
      const a = $("dFrom").value;
      const start = a ? tms(`${a} 00:00`) : tms(state.from);
      state.from = tstr(start);
      state.to = tstr(start + +b.dataset.days * DAY);
      render();
    });
    $("resSel").onchange = (e) => { state.res = e.target.value; render(); };
    document.querySelectorAll(".seg button").forEach((b) => b.onclick = () => { state.mode = b.dataset.mode; render(); });
    $("logSuction").onchange = (e) => { state.logSuction = e.target.checked; render(); };
    $("splitRain").onchange = (e) => { state.splitRain = e.target.checked; render(); };
    $("exportBtn").onclick = exportCSV;
    $("themeBtn").onclick = () => {
      const dark = getComputedStyle(document.documentElement).colorScheme === "dark";
      document.documentElement.dataset.theme = dark ? "light" : "dark";
      renderSeriesList(); render();
    };
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { renderSeriesList(); render(); });
  }

  // ---------------------------------------------------------------- avvio
  function registerItalian() {
    Plotly.register({
      moduleType: "locale", name: "it", dictionary: {},
      format: {
        days: ["Domenica", "Lunedì", "Martedì", "Mercoledì", "Giovedì", "Venerdì", "Sabato"],
        shortDays: ["Dom", "Lun", "Mar", "Mer", "Gio", "Ven", "Sab"],
        months: ["Gennaio", "Febbraio", "Marzo", "Aprile", "Maggio", "Giugno", "Luglio", "Agosto",
                 "Settembre", "Ottobre", "Novembre", "Dicembre"],
        shortMonths: ["Gen", "Feb", "Mar", "Apr", "Mag", "Giu", "Lug", "Ago", "Set", "Ott", "Nov", "Dic"],
        date: "%d/%m/%Y", decimal: ".", thousands: " ",
      },
    });
  }
  async function init() {
    if (window.Plotly) registerItalian();
    try {
      catalog = await fetch("catalog.json").then((r) => r.json());
    } catch (e) {
      $("chart").innerHTML = `<div class="empty">Impossibile leggere catalog.json (${e.message}).<br>
        Il portale va aperto da un server web (es. GitHub Pages), non come file locale.</div>`;
      return;
    }
    catalog.series.forEach((s) => (seriesById[s.series_id] = s));
    catalog.stations.forEach((s) => (stationById[s.station_id] = s));
    readHash();
    bindControls();
    renderSeriesList();
    initMap();
    render();
  }
  init();
})();
