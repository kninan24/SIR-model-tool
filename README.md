# SIR Model Tool

Static GitHub Pages dashboard for an educational influenza SIR model.

The app runs fully in the browser using HTML, CSS, JavaScript, and Chart.js. It compares a simulated SIR model against an embedded hypothetical weekly influenza-like illness dataset and lets users adjust model parameters with sliders.

Repository target: <https://github.com/kninan24/SIR-model-tool>

## Files

```text
index.html
assets/styles.css
assets/app.js
.nojekyll
README.md
```

## Model

The dashboard simulates a deterministic SIR model:

```text
dS/dt = -beta * S * I / N
dI/dt = beta * S * I / N - gamma * I
dR/dt = gamma * I
```

Parameters:

- `N`: population size
- `I0`: initial infectious individuals
- `R0`: initially immune/recovered individuals
- `beta`: transmission rate per day
- `gamma`: recovery rate per day
- `R0_basic`: basic reproduction number, where `beta = R0_basic * gamma`
- `rho`: reporting/scaling multiplier used to compare weekly SIR incidence against the embedded observed data scale

The numerical simulation uses a fourth-order Runge-Kutta method in JavaScript.

## Dashboard features

- SIR model sliders for population, initial infections, initial immunity, R0, infectious period, simulation duration, and reporting multiplier
- Embedded hypothetical weekly observed influenza-like illness data
- Model-vs-data comparison plot
- S/I/R compartment plot
- Daily incidence and effective Rt plot
- RMSE metric showing how close the current slider settings are to the hypothetical observed data
- CSV download of weekly model-vs-data comparison

## Local preview

Because this is a static site, you can open `index.html` directly in a browser.

For a cleaner local preview, from this folder run:

```bash
python -m http.server 8000
```

Then open:

```text
http://localhost:8000
```

## GitHub repository setup

From the parent folder, run:

```bash
cd SIR-model-tool
git status
git add index.html assets/styles.css assets/app.js README.md .nojekyll
git commit -m "Create static SIR model dashboard tool"
git remote add origin https://github.com/kninan24/SIR-model-tool.git
git branch -M main
git push -u origin main
```

If the remote already exists, use:

```bash
git remote set-url origin https://github.com/kninan24/SIR-model-tool.git
git push -u origin main
```

## GitHub Pages setup

On GitHub:

1. Open the repository.
2. Go to **Settings**.
3. Go to **Pages**.
4. Under **Build and deployment**, choose:
   - Source: **Deploy from a branch**
   - Branch: **main**
   - Folder: `/root` if this folder is the repository root.
5. Click **Save**.

For this repository, use the repository root as the GitHub Pages folder because `index.html` is at the root of `SIR-model-tool`.

## Notes

This is an educational dashboard. It should not be interpreted as a real influenza forecast.
