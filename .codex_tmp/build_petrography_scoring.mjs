import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const repo = "/Users/allen/Documents/GitHub/bighorn_paleogene_proxies";
const imageRoot = "/Users/allen/Library/CloudStorage/OneDrive-Personal/MLA Work/NSF_EARPF/Images/Thin Sections";
const outputDir = path.join(repo, "outputs", "petrography_scoring");

function parseCsv(text) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (ch === '"') quoted = false;
      else field += ch;
    } else {
      if (ch === '"') quoted = true;
      else if (ch === ',') { row.push(field); field = ""; }
      else if (ch === '\n') { row.push(field.replace(/\r$/, "")); rows.push(row); row = []; field = ""; }
      else field += ch;
    }
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  const headers = rows.shift();
  return rows.filter(r => r.some(x => x !== "")).map(r => Object.fromEntries(headers.map((h, i) => [h, r[i] ?? ""])));
}

async function walk(dir) {
  const out = [];
  for (const ent of await fs.readdir(dir, { withFileTypes: true })) {
    const p = path.join(dir, ent.name);
    if (ent.isDirectory()) out.push(...await walk(p)); else out.push(p);
  }
  return out;
}

const sampleDirs = (await fs.readdir(imageRoot, { withFileTypes: true }))
  .filter(x => x.isDirectory() && x.name.startsWith("PK95"))
  .map(x => x.name)
  .sort();

const summaries = [];
const empaRows = [];
for (const sample of sampleDirs) {
  const folder = path.join(imageRoot, sample);
  const files = await walk(folder);
  const lower = files.map(f => path.basename(f).toLowerCase());
  const aois = [...new Set(files.map(f => {
    const m = f.match(/(PK95[^/]*_AOI\d+)/i); return m ? m[1] : null;
  }).filter(Boolean))].sort();
  const count = re => lower.filter(n => re.test(n)).length;
  summaries.push({
    sample,
    folder,
    ppl: count(/ppl/),
    xpl: count(/xpl/),
    marked: count(/marked/),
    bse: lower.filter(n => n.includes("bse") && n.endsWith(".png")).length,
    cl: lower.filter(n => /_cl/.test(n) && n.endsWith(".png")).length,
    maps: ["mn", "fe", "mg", "si", "ca"].filter(e => lower.some(n => n.includes(`_${e}`) && n.endsWith(".png"))).join(", "),
    aois: aois.length,
  });
  for (const aoi of aois) empaRows.push([sample, aoi]);
}

const csv = await fs.readFile(path.join(repo, "data/processed/CFB_soilwater_reconstruction_summary_age_calibrated.csv"), "utf8");
const isotope = parseCsv(csv);
const normalize = s => s === "PK95-SC-80.1" ? "PK95-SC-80" : s;
const isotopeById = new Map(isotope.map(r => [r.MLA_horizon_id, r]));

const wb = Workbook.create();
const instructions = wb.worksheets.add("Read me");
const petro = wb.worksheets.add("Petrography scoring");
const empa = wb.worksheets.add("EMPA scoring");
const synthesis = wb.worksheets.add("Synthesis");
const context = wb.worksheets.add("Isotope context");

const navy = "#183B56", teal = "#2A7F86", pale = "#EAF3F4", gold = "#D9A441", light = "#F5F7F8", red = "#F4CCCC";
for (const sh of [instructions, petro, empa, synthesis, context]) sh.showGridLines = false;

instructions.getRange("A1:H1").merge();
instructions.getRange("A1").values = [["Bighorn Basin carbonate preservation scoring"]];
instructions.getRange("A1:H1").format = { fill: navy, font: { color: "#FFFFFF", bold: true, size: 18 }, rowHeight: 32, verticalAlignment: "center" };
instructions.getRange("A3:H3").merge();
instructions.getRange("A3").values = [["Purpose"]];
instructions.getRange("A4:H5").merge();
instructions.getRange("A4").values = [["Record observations first, then infer late-alteration risk. Crystal size and recrystallization texture do not establish timing: micrite, microspar, and spar may each be early or late. Estimate abundance only for carbonate shown to be late by crosscutting, replacement, CL, or elemental relationships. Petrographic and EMPA AOIs are not spatially identical and must be integrated only at sample level."]];
instructions.getRange("A7:B13").values = [
  ["Evidence category", "Definition"],
  ["Matrix texture", "Descriptive grain-size/fabric observation only; never use micrite versus microspar alone to assign timing."],
  ["Demonstrably late carbonate", "Carbonate with clear crosscutting, pore/fracture-filling, truncation, or replacement relationships."],
  ["Timing unresolved", "A texture or chemical domain is present, but its age relative to pedogenic carbonate cannot be established."],
  ["Localized late alteration", "Clearly late carbonate is spatially restricted; record its observed extent without extrapolating beyond the imaged area."],
  ["Pervasive late alteration", "Clearly late replacement or cement is widespread and fabric-destructive across the examined material."],
  ["Abundance bins", "Use <10%, 10–30%, 30–60%, or >60% only for demonstrably late carbonate; otherwise choose not estimable."],
];
instructions.getRange("A14:H14").merge();
instructions.getRange("A14").values = [["Workflow: record petrography without consulting isotope values → assess EMPA independently → combine timing evidence in Synthesis → consult Isotope context only after alteration-risk decisions are locked."]];
instructions.getRange("A3:H3").format = { fill: teal, font: { color: "#FFFFFF", bold: true } };
instructions.getRange("A7:B7").format = { fill: teal, font: { color: "#FFFFFF", bold: true } };
instructions.getRange("A14:H14").format = { fill: "#FFF2CC", font: { bold: true }, wrapText: true, rowHeight: 38 };
instructions.getRange("A4:H5").format.wrapText = true;
instructions.getRange("A7:B13").format.wrapText = true;
instructions.getRange("A1:H14").format.font = { name: "Arial", size: 10 };
instructions.getRange("A1:H1").format.font = { name: "Arial", size: 18, bold: true, color: "#FFFFFF" };
instructions.getRange("A:A").format.columnWidth = 22;
instructions.getRange("B:B").format.columnWidth = 82;

const petroHeaders = ["Sample ID","Folder","PPL fields","XPL fields","Marked image","Dominant matrix texture","Organized/pedogenic fabric extent","Fabric description","Clearly crosscutting cement extent","Replacement/fabric destruction evidence","Veins/fractures extent","Demonstrably late carbonate abundance","Petrographic timing inference","Coverage adequate?","Confidence","Microscope review complete","Notes"];
petro.getRange("A1:Q1").values = [petroHeaders];
petro.getRange(`A2:Q${summaries.length + 1}`).values = summaries.map(s => [s.sample,s.folder,s.ppl,s.xpl,s.marked > 0 ? "Yes" : "No","","","","","","","","","","","No",""]);
petro.getRange("A1:Q1").format = { fill: navy, font: { color: "#FFFFFF", bold: true }, wrapText: true, rowHeight: 42 };
petro.getRange(`A2:Q${summaries.length + 1}`).format = { font: { name: "Arial", size: 9 }, verticalAlignment: "top" };
petro.freezePanes.freezeRows(1); petro.freezePanes.freezeColumns(1);
petro.tables.add(`A1:Q${summaries.length + 1}`, true, "PetrographyScores").style = "TableStyleMedium2";
petro.getRange(`F2:F${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["micrite","microspar","mixed micrite–microspar","coarse crystalline","indeterminate"] } };
petro.getRange(`G2:G${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["absent","localized","widespread","indeterminate"] } };
petro.getRange(`I2:I${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["absent","localized","distributed","pervasive","indeterminate"] } };
petro.getRange(`J2:J${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["absent","suspected","clear localized","clear pervasive","indeterminate"] } };
petro.getRange(`K2:K${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["absent","minor","moderate","abundant","indeterminate"] } };
petro.getRange(`L2:L${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["<10%","10–30%","30–60%",">60%","not estimable"] } };
petro.getRange(`M2:M${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["no clear petrographic evidence of late matrix alteration","localized demonstrably late carbonate","pervasive demonstrably late alteration","timing unresolved"] } };
petro.getRange(`N2:N${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["yes","no","uncertain"] } };
petro.getRange(`O2:O${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["high","moderate","low"] } };
petro.getRange(`P2:P${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["Yes","No"] } };
petro.getRange(`M2:M${summaries.length + 1}`).conditionalFormats.add("containsText", { text: "no clear", format: { fill: "#D9EAD3" } });
petro.getRange(`M2:M${summaries.length + 1}`).conditionalFormats.add("containsText", { text: "localized", format: { fill: "#FFF2CC" } });
petro.getRange(`M2:M${summaries.length + 1}`).conditionalFormats.add("containsText", { text: "pervasive", format: { fill: red } });

const empaHeaders = ["Sample ID","EMPA AOI","BSE texture/heterogeneity","CL pattern","Elemental-domain extent","Domain relationship to fractures or fabrics","Evidence domain is late?","EMPA timing inference","AOI coverage adequate?","Confidence","Notes"];
empa.getRange("A1:K1").values = [empaHeaders];
empa.getRange(`A2:K${empaRows.length + 1}`).values = empaRows.map(r => [...r,"","","","","","","","",""]);
empa.getRange("A1:K1").format = { fill: teal, font: { color: "#FFFFFF", bold: true }, wrapText: true, rowHeight: 42 };
empa.freezePanes.freezeRows(1); empa.freezePanes.freezeColumns(2);
empa.tables.add(`A1:K${empaRows.length + 1}`, true, "EMPAScores").style = "TableStyleMedium4";
empa.getRange(`C2:C${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["homogeneous","localized heterogeneity","distributed heterogeneity","indeterminate"] } };
empa.getRange(`D2:D${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["uniform","patchy","zoned","fracture-associated","indeterminate"] } };
empa.getRange(`E2:E${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["absent","localized","distributed","pervasive","indeterminate"] } };
empa.getRange(`G2:G${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["yes","no","uncertain"] } };
empa.getRange(`H2:H${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["no clear evidence of late alteration","localized late chemical domain","distributed/pervasive late chemical alteration","chemical heterogeneity; timing unresolved"] } };
empa.getRange(`I2:I${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["yes","no","uncertain"] } };
empa.getRange(`J2:J${empaRows.length + 1}`).dataValidation = { rule: { type: "list", values: ["high","moderate","low"] } };

const synthHeaders = ["Sample ID","Petrographic timing inference","Demonstrably late carbonate abundance","EMPA timing inference","AOI relationship","Combined late-alteration risk","Primary-record decision","Decision basis","Locked before isotope review?"];
synthesis.getRange("A1:I1").values = [synthHeaders];
synthesis.getRange(`A2:I${summaries.length + 1}`).values = summaries.map(s => [s.sample,"","","","Different AOIs—sample-level integration only","","","",""]);
synthesis.getRange("A1:J1").format = { fill: gold, font: { color: "#FFFFFF", bold: true }, wrapText: true, rowHeight: 42 };
synthesis.freezePanes.freezeRows(1); synthesis.freezePanes.freezeColumns(1);
synthesis.tables.add(`A1:I${summaries.length + 1}`, true, "PreservationSynthesis").style = "TableStyleMedium9";
synthesis.getRange(`B2:B${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["no clear petrographic evidence of late matrix alteration","localized demonstrably late carbonate","pervasive demonstrably late alteration","timing unresolved"] } };
synthesis.getRange(`C2:C${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["<10%","10–30%","30–60%",">60%","not estimable"] } };
synthesis.getRange(`D2:D${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["no clear evidence of late alteration","localized late chemical domain","distributed/pervasive late chemical alteration","chemical heterogeneity; timing unresolved","No EMPA"] } };
synthesis.getRange(`F2:F${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["low observed risk","localized alteration risk","high/pervasive alteration risk","unresolved"] } };
synthesis.getRange(`G2:G${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["Retain","Sensitivity-only","Exclude","Pending"] } };
synthesis.getRange(`I2:I${summaries.length + 1}`).dataValidation = { rule: { type: "list", values: ["Yes","No"] } };

const ctxHeaders = ["Sample ID","Strat height (m)","Age (Ma)","U-M T47 (°C)","U-M T47 SE","Carbonate δ18O (VSMOW)","Carbonate Δ′17O (per meg)","Δ′17O primary use","Current T screen exclusion","Geochemical preservation index (not petrographic probability)","Image coverage","EMPA AOIs"];
context.getRange("A1:L1").values = [ctxHeaders];
context.getRange(`A2:L${summaries.length + 1}`).values = summaries.map(s => {
  const id = normalize(s.sample); const r = isotopeById.get(id) ?? {};
  const val = k => (r[k] && r[k] !== "NA" ? (Number.isFinite(Number(r[k])) ? Number(r[k]) : r[k]) : null);
  return [s.sample,val("strat_height_m"),val("Age_Ma"),val("IPLD47_mean_T47_C"),val("IPLD47_se_T47_C"),val("IPL_NuDog_d18Ocarb_VSMOW"),val("IPL17O_mean_Dp17Ocarb"),val("D17O_primary_use"),val("exclude_from_temp_model"),val("p_altered_preservation.x"),`PPL ${s.ppl}; XPL ${s.xpl}; marked ${s.marked > 0 ? "yes" : "no"}`,s.aois];
});
context.getRange("A1:L1").format = { fill: "#666666", font: { color: "#FFFFFF", bold: true }, wrapText: true, rowHeight: 40 };
context.freezePanes.freezeRows(1); context.freezePanes.freezeColumns(1);
context.tables.add(`A1:L${summaries.length + 1}`, true, "IsotopeContext").style = "TableStyleMedium3";
context.getRange(`B2:F${summaries.length + 1}`).format.numberFormat = "0.00";
context.getRange(`G2:G${summaries.length + 1}`).format.numberFormat = "0.0";
context.getRange(`J2:J${summaries.length + 1}`).format.numberFormat = "0.00";

for (const sh of [petro, empa, synthesis, context]) {
  const used = sh.getUsedRange(); used.format.font = { name: "Arial", size: 9 };
  sh.getRange("1:1").format.font = { name: "Arial", size: 9, bold: true, color: "#FFFFFF" };
  used.format.autofitColumns();
  used.format.autofitRows();
}
petro.getRange("A:A").format.columnWidth = 18; petro.getRange("B:B").format.columnWidth = 42; petro.getRange("F:P").format.columnWidth = 22; petro.getRange("Q:Q").format.columnWidth = 42;
empa.getRange("A:B").format.columnWidth = 20; empa.getRange("C:J").format.columnWidth = 24; empa.getRange("K:K").format.columnWidth = 42;
synthesis.getRange("A:A").format.columnWidth = 18; synthesis.getRange("B:I").format.columnWidth = 25;
context.getRange("A:A").format.columnWidth = 18; context.getRange("B:J").format.columnWidth = 18; context.getRange("K:K").format.columnWidth = 30;

await fs.mkdir(outputDir, { recursive: true });
const outputPath = path.join(outputDir, "BHB_petrography_EMPA_scoring_REVISED.xlsx");
const file = await SpreadsheetFile.exportXlsx(wb);
await file.save(outputPath);

for (const sh of [instructions, petro, empa, synthesis, context]) {
  const preview = await wb.render({ sheetName: sh.name, autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(path.join(outputDir, `preview_${sh.name.replaceAll(" ", "_")}.png`), new Uint8Array(await preview.arrayBuffer()));
}

const check = await wb.inspect({ kind: "table", sheetId: "Petrography scoring", range: `A1:Q${summaries.length + 1}`, include: "values,formulas", tableMaxRows: 8, tableMaxCols: 17, maxChars: 7000 });
const errors = await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
console.log(check.ndjson);
console.log(errors.ndjson);
console.log(outputPath);
