const observedWeeklyData = [
  { week: 0, ili: 0.8 }, { week: 1, ili: 0.9 }, { week: 2, ili: 1.1 },
  { week: 3, ili: 1.4 }, { week: 4, ili: 2.0 }, { week: 5, ili: 3.1 },
  { week: 6, ili: 4.8 }, { week: 7, ili: 6.4 }, { week: 8, ili: 7.2 },
  { week: 9, ili: 6.5 }, { week: 10, ili: 5.1 }, { week: 11, ili: 3.7 },
  { week: 12, ili: 2.6 }, { week: 13, ili: 1.9 }, { week: 14, ili: 1.4 },
  { week: 15, ili: 1.1 }, { week: 16, ili: 0.9 }, { week: 17, ili: 0.8 },
  { week: 18, ili: 0.7 }, { week: 19, ili: 0.7 }, { week: 20, ili: 0.6 },
  { week: 21, ili: 0.6 }, { week: 22, ili: 0.6 }, { week: 23, ili: 0.7 },
  { week: 24, ili: 0.8 }, { week: 25, ili: 0.9 }
];

const presets = {
  mild: { population: 10000, initialInfectious: 15, initialImmune: 15, r0: 1.20, infectiousPeriod: 4.0, duration: 180, reportingScale: 0.020 },
  typical: { population: 10000, initialInfectious: 25, initialImmune: 10, r0: 1.45, infectiousPeriod: 4.0, duration: 180, reportingScale: 0.020 },
  severe: { population: 10000, initialInfectious: 40, initialImmune: 5, r0: 1.90, infectiousPeriod: 3.5, duration: 180, reportingScale: 0.012 }
};

const controls = {
  preset: document.getElementById('preset'),
  population: document.getElementById('population'),
  initialInfectious: document.getElementById('initialInfectious'),
  initialImmune: document.getElementById('initialImmune'),
  r0: document.getElementById('r0'),
  infectiousPeriod: document.getElementById('infectiousPeriod'),
  duration: document.getElementById('duration'),
  reportingScale: document.getElementById('reportingScale')
};

const outputs = {
  population: document.getElementById('populationValue'),
  initialInfectious: document.getElementById('initialInfectiousValue'),
  initialImmune: document.getElementById('initialImmuneValue'),
  r0: document.getElementById('r0Value'),
  infectiousPeriod: document.getElementById('infectiousPeriodValue'),
  duration: document.getElementById('durationValue'),
  reportingScale: document.getElementById('reportingScaleValue'),
  beta: document.getElementById('betaMetric'),
  gamma: document.getElementById('gammaMetric'),
  peak: document.getElementById('peakMetric'),
  attack: document.getElementById('attackMetric'),
  rmse: document.getElementById('rmseMetric')
};

const resetButton = document.getElementById('resetButton');
const downloadButton = document.getElementById('downloadButton');
let fitChart, compartmentChart, incidenceChart;
let latestWeeklyRows = [];

function numberFormat(value, digits = 0) {
  return Number(value).toLocaleString(undefined, { maximumFractionDigits: digits, minimumFractionDigits: digits });
}

function getParams() {
  const N = Number(controls.population.value);
  const I0 = Number(controls.initialInfectious.value);
  const immunePct = Number(controls.initialImmune.value) / 100;
  const R0Initial = Math.round(N * immunePct);
  const S0 = Math.max(0, N - I0 - R0Initial);
  const r0 = Number(controls.r0.value);
  const infectiousPeriod = Number(controls.infectiousPeriod.value);
  const gamma = 1 / infectiousPeriod;
  const beta = r0 * gamma;
  const duration = Number(controls.duration.value);
  const reportingScale = Number(controls.reportingScale.value);
  return { N, S0, I0, R0Initial, r0, infectiousPeriod, gamma, beta, duration, reportingScale };
}

function updateControlLabels(p) {
  outputs.population.textContent = numberFormat(p.N);
  outputs.initialInfectious.textContent = numberFormat(p.I0);
  outputs.initialImmune.textContent = `${numberFormat(Number(controls.initialImmune.value), 0)}%`;
  outputs.r0.textContent = p.r0.toFixed(2);
  outputs.infectiousPeriod.textContent = p.infectiousPeriod.toFixed(1);
  outputs.duration.textContent = numberFormat(p.duration);
  outputs.reportingScale.textContent = p.reportingScale.toFixed(3);
}

function sirDerivatives(state, p) {
  const incidence = p.beta * state.S * state.I / p.N;
  const recoveries = p.gamma * state.I;
  return { dS: -incidence, dI: incidence - recoveries, dR: recoveries, incidence };
}

function rk4Step(state, p, dt) {
  const k1 = sirDerivatives(state, p);
  const k2 = sirDerivatives({ S: state.S + 0.5 * dt * k1.dS, I: state.I + 0.5 * dt * k1.dI, R: state.R + 0.5 * dt * k1.dR }, p);
  const k3 = sirDerivatives({ S: state.S + 0.5 * dt * k2.dS, I: state.I + 0.5 * dt * k2.dI, R: state.R + 0.5 * dt * k2.dR }, p);
  const k4 = sirDerivatives({ S: state.S + dt * k3.dS, I: state.I + dt * k3.dI, R: state.R + dt * k3.dR }, p);
  return {
    S: Math.max(0, state.S + (dt / 6) * (k1.dS + 2 * k2.dS + 2 * k3.dS + k4.dS)),
    I: Math.max(0, state.I + (dt / 6) * (k1.dI + 2 * k2.dI + 2 * k3.dI + k4.dI)),
    R: Math.max(0, state.R + (dt / 6) * (k1.dR + 2 * k2.dR + 2 * k3.dR + k4.dR))
  };
}

function simulateSIR(p) {
  const dt = 0.25;
  const rows = [];
  let state = { S: p.S0, I: p.I0, R: p.R0Initial };
  for (let day = 0; day <= p.duration; day += 1) {
    const deriv = sirDerivatives(state, p);
    rows.push({ day, S: state.S, I: state.I, R: state.R, incidence: deriv.incidence, Rt: p.r0 * (state.S / p.N) });
    for (let sub = 0; sub < 1 / dt; sub += 1) state = rk4Step(state, p, dt);
  }
  return rows;
}

function weeklyComparison(rows, p) {
  return observedWeeklyData.map(obs => {
    const start = obs.week * 7;
    const end = start + 6;
    const weekRows = rows.filter(row => row.day >= start && row.day <= end);
    const modelIncidence = weekRows.reduce((sum, row) => sum + row.incidence, 0);
    const modelILI = p.reportingScale * modelIncidence;
    return { week: obs.week, observed_ili: obs.ili, model_ili: modelILI, weekly_incidence: modelIncidence, residual: obs.ili - modelILI };
  });
}

function rmse(rows) {
  const mse = rows.reduce((sum, row) => sum + row.residual ** 2, 0) / rows.length;
  return Math.sqrt(mse);
}

function summarize(rows, p) {
  const peak = rows.reduce((best, row) => row.I > best.I ? row : best, rows[0]);
  const finalAttackRate = ((p.S0 - rows[rows.length - 1].S) / p.N) * 100;
  return { peak, finalAttackRate };
}

function baseOptions(xTitle, yTitle) {
  return {
    responsive: true,
    maintainAspectRatio: false,
    interaction: { mode: 'index', intersect: false },
    plugins: { legend: { position: 'top' }, tooltip: { callbacks: { label: ctx => `${ctx.dataset.label}: ${numberFormat(ctx.parsed.y, 2)}` } } },
    scales: { x: { title: { display: true, text: xTitle }, grid: { display: false } }, y: { title: { display: true, text: yTitle }, beginAtZero: true } }
  };
}

function renderCharts(rows, weeklyRows) {
  const dayLabels = rows.map(row => row.day);
  const weekLabels = weeklyRows.map(row => `Week ${row.week}`);
  if (fitChart) fitChart.destroy();
  if (compartmentChart) compartmentChart.destroy();
  if (incidenceChart) incidenceChart.destroy();

  fitChart = new Chart(document.getElementById('fitChart'), {
    type: 'line',
    data: { labels: weekLabels, datasets: [
      { label: 'Hypothetical observed ILI', data: weeklyRows.map(row => row.observed_ili), borderColor: '#0072B2', backgroundColor: '#0072B222', borderWidth: 2.5, pointRadius: 4 },
      { label: 'Scaled SIR weekly incidence', data: weeklyRows.map(row => row.model_ili), borderColor: '#D55E00', backgroundColor: '#D55E0022', borderWidth: 2.5, pointRadius: 3 }
    ] },
    options: baseOptions('Epidemic week', 'ILI-like value')
  });

  compartmentChart = new Chart(document.getElementById('compartmentChart'), {
    type: 'line',
    data: { labels: dayLabels, datasets: [
      { label: 'Susceptible', data: rows.map(row => row.S), borderColor: '#0072B2', backgroundColor: '#0072B222', borderWidth: 2.5, pointRadius: 0 },
      { label: 'Infectious', data: rows.map(row => row.I), borderColor: '#D55E00', backgroundColor: '#D55E0022', borderWidth: 2.5, pointRadius: 0 },
      { label: 'Recovered', data: rows.map(row => row.R), borderColor: '#009E73', backgroundColor: '#009E7322', borderWidth: 2.5, pointRadius: 0 }
    ] },
    options: baseOptions('Day', 'Number of individuals')
  });

  incidenceChart = new Chart(document.getElementById('incidenceChart'), {
    type: 'line',
    data: { labels: dayLabels, datasets: [
      { label: 'Daily incidence', data: rows.map(row => row.incidence), borderColor: '#CC79A7', backgroundColor: '#CC79A722', borderWidth: 2.5, pointRadius: 0, yAxisID: 'y' },
      { label: 'Effective Rt', data: rows.map(row => row.Rt), borderColor: '#7B2CBF', backgroundColor: '#7B2CBF22', borderWidth: 2.5, pointRadius: 0, yAxisID: 'y1' }
    ] },
    options: { ...baseOptions('Day', 'New infections per day'), scales: { x: { title: { display: true, text: 'Day' }, grid: { display: false } }, y: { title: { display: true, text: 'New infections per day' }, beginAtZero: true, position: 'left' }, y1: { title: { display: true, text: 'Effective Rt' }, beginAtZero: true, position: 'right', grid: { drawOnChartArea: false } } } }
  });
}

function updateDashboard() {
  let p = getParams();
  if (p.I0 + p.R0Initial >= p.N) controls.initialInfectious.value = Math.max(1, p.N - p.R0Initial - 1);
  p = getParams();
  updateControlLabels(p);

  const rows = simulateSIR(p);
  latestWeeklyRows = weeklyComparison(rows, p);
  const summary = summarize(rows, p);
  const fitRmse = rmse(latestWeeklyRows);

  outputs.beta.textContent = `${p.beta.toFixed(3)}/day`;
  outputs.gamma.textContent = `${p.gamma.toFixed(3)}/day`;
  outputs.peak.textContent = `${numberFormat(summary.peak.I, 0)} on day ${summary.peak.day}`;
  outputs.attack.textContent = `${summary.finalAttackRate.toFixed(1)}%`;
  outputs.rmse.textContent = fitRmse.toFixed(2);

  renderCharts(rows, latestWeeklyRows);
}

function applyPreset(name) {
  const preset = presets[name];
  if (!preset) return;
  controls.population.value = preset.population;
  controls.initialInfectious.value = preset.initialInfectious;
  controls.initialImmune.value = preset.initialImmune;
  controls.r0.value = preset.r0;
  controls.infectiousPeriod.value = preset.infectiousPeriod;
  controls.duration.value = preset.duration;
  controls.reportingScale.value = preset.reportingScale;
  controls.preset.value = name;
  updateDashboard();
}

function downloadCSV() {
  if (!latestWeeklyRows.length) return;
  const header = Object.keys(latestWeeklyRows[0]);
  const csv = [
    header.join(','),
    ...latestWeeklyRows.map(row => header.map(key => Number(row[key]).toFixed(6)).join(','))
  ].join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = 'sir_model_vs_hypothetical_data.csv';
  link.click();
  URL.revokeObjectURL(url);
}

Object.values(controls).forEach(control => {
  control.addEventListener('input', () => {
    if (control !== controls.preset) controls.preset.value = 'custom';
    updateDashboard();
  });
});

controls.preset.addEventListener('change', event => {
  if (event.target.value !== 'custom') applyPreset(event.target.value);
});

resetButton.addEventListener('click', () => applyPreset('typical'));
downloadButton.addEventListener('click', downloadCSV);
window.addEventListener('DOMContentLoaded', () => applyPreset('typical'));
