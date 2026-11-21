# Visualizador de la estimación de la KAN en las estaciones de prueba (1998-2003)
# Uso: uv run streamlit run visualizador.py
import calendar
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly
import streamlit as st
import streamlit.components.v1 as components

RAIZ = Path(__file__).parent
MESES = ['ene', 'feb', 'mar', 'abr', 'may', 'jun', 'jul', 'ago', 'sep', 'oct', 'nov', 'dic']


# Fase ENSO del año con el ONI promedio (umbral de 0,5, criterio NOAA)
def fase_enso(oni):
    return 'El Niño' if oni >= 0.5 else ('La Niña' if oni <= -0.5 else 'Neutral')


# Días de cada péntada en un año: 6 por mes, la última va hasta fin de mes
def etiquetas_pentadas(anio):
    etiquetas = []
    for k in range(72):
        mes, j = k // 6 + 1, k % 6
        inicio = 1 + 5 * j
        fin = calendar.monthrange(anio, mes)[1] if j == 5 else inicio + 4
        etiquetas.append(f'Péntada {k + 1} · {inicio}–{fin} {MESES[mes - 1]}')
    return etiquetas


@st.cache_data
def cargar_datos():
    meta = pd.read_csv(RAIZ / 'Curvas/metadatos_curvas.csv', dtype={'curva': str, 'STATION_ID': str})
    leer = lambda f: pd.read_csv(RAIZ / 'Curvas' / f, dtype={'curva': str}).set_index('curva')
    estaciones, satelite = leer('curvas_estaciones_72_pentadas.csv').clip(lower=0), leer('curvas_satelite_72_pentadas.csv').clip(lower=0)
    kan = pd.read_csv(RAIZ / 'predicciones_estaciones_reservadas.csv', dtype={'station_id': str})
    oni = pd.read_csv(RAIZ / 'data/oni_anual.csv')

    # Curvas de prueba: KAN, curva suavizada de la estación y satélite, con la RECM de la KAN
    prueba = meta[meta.particion == 'prueba']
    por_anio = {}
    for anio, grupo in prueba.groupby('YEAR'):
        lista = []
        for _, fila in grupo.iterrows():
            pred = kan[(kan.station_id == fila.STATION_ID) & (kan.year == anio)].sort_values('day').value_mm.to_numpy()
            est = estaciones.loc[fila.curva].to_numpy()
            lista.append({'id': fila.STATION_ID, 'lon': fila.longitud, 'lat': fila.latitud,
                          'altitud': round(float(fila.altitud_completa)), 'recm': round(float(np.sqrt(np.mean((pred - est) ** 2))), 1),
                          'kan': np.round(pred, 1).tolist(), 'estacion': np.round(est, 1).tolist(),
                          'satelite': np.round(satelite.loc[fila.curva].to_numpy(), 1).tolist()})
        por_anio[int(anio)] = {'estaciones': lista, 'pentadas': etiquetas_pentadas(int(anio))}

    # Estaciones de entrenamiento, solo como referencia en el mapa
    contexto = meta[meta.particion == 'entrenamiento'].drop_duplicates('STATION_ID')[['longitud', 'latitud']]
    return {'por_anio': por_anio,
            'contexto': {'lon': contexto.longitud.round(4).tolist(), 'lat': contexto.latitud.round(4).tolist()},
            'oni': {'anio': oni.YEAR.tolist(), 'valor': oni.oni_promedio.round(2).tolist(),
                    'fase': [fase_enso(v) for v in oni.oni_promedio]}}


st.set_page_config(page_title='Estimación de la KAN', layout='wide')
st.markdown('#### Estimación de la lluvia por péntada en las estaciones de prueba')
st.caption('Arrastra la barra de años o haz clic en un año de prueba del ONI. Pasa el cursor por una estación para ver su curva '
           'y por la curva para ver la péntada. Las estaciones de prueba no se usaron para entrenar el modelo.')
relieve = st.toggle('Relieve', value=True)

datos = cargar_datos()

HTML = """
<style>
  body { font-family: sans-serif; margin: 0; color: #2c2c2a; }
  .fila { display: flex; gap: 12px; }
  .caja { border: 0.5px solid #d3d1c7; border-radius: 12px; padding: 8px; }
  .barra { display: flex; align-items: center; gap: 10px; font-size: 13px; color: #5f5e5a; margin: 2px 6px 0; }
  input[type=range] { flex: 1; }
</style>
<div class="caja">
  <div id="oni" style="height:130px"></div>
  <div class="barra"><span>1998</span><input id="anio" type="range" min="1998" max="2003" step="1" value="1998"><span>2003</span>
    <b id="etiqueta" style="min-width:170px; text-align:right"></b></div>
</div>
<div class="fila" style="margin-top:12px">
  <div class="caja" style="flex:1"><div id="mapa" style="height:520px"></div></div>
  <div class="caja" style="flex:1"><div id="curva" style="height:520px"></div></div>
</div>
<script>PLOTLYJS</script>
<script>
const D = DATOS;
const RELIEVE = RELIEVE_ACTIVO;
const fuente = {family: 'sans-serif', size: 12, color: '#2c2c2a'};
let anio = 1998, actual = null;

function faseAnio(a) {
  const i = D.oni.anio.indexOf(a);
  return D.oni.fase[i] + ' · ONI ' + D.oni.valor[i].toFixed(2).replace('.', ',');
}

// Serie del ONI promedio; los años de prueba se marcan y se pueden elegir con un clic
function dibujarOni() {
  const prueba = D.oni.anio.map(a => a >= 1998);
  Plotly.react('oni', [
    {x: D.oni.anio, y: D.oni.valor, type: 'scatter', mode: 'lines', line: {color: '#888780', width: 1.5}, hoverinfo: 'skip'},
    {x: D.oni.anio.filter((a, i) => prueba[i]), y: D.oni.valor.filter((v, i) => prueba[i]), type: 'scatter', mode: 'markers',
     marker: {size: 9, color: D.oni.anio.filter(a => a >= 1998).map(a => a === anio ? '#E24B4A' : '#B4B2A9')},
     hovertemplate: '%{x}: ONI %{y:.2f}<extra></extra>'}
  ], {
    margin: {l: 40, r: 10, t: 22, b: 24}, showlegend: false, font: fuente, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    title: {text: 'ONI promedio anual', font: {size: 12}, x: 0.01, y: 0.97},
    xaxis: {dtick: 2, showgrid: false}, yaxis: {zeroline: true, zerolinecolor: '#d3d1c7', gridcolor: '#f1efe8'},
    shapes: [
      {type: 'rect', x0: 1997.5, x1: 2003.5, y0: 0, y1: 1, yref: 'paper', fillcolor: '#f1efe8', line: {width: 0}, layer: 'below'},
      {type: 'line', x0: 1980.5, x1: 2003.5, y0: 0.5, y1: 0.5, line: {color: '#E24B4A', width: 0.8, dash: 'dot'}},
      {type: 'line', x0: 1980.5, x1: 2003.5, y0: -0.5, y1: -0.5, line: {color: '#378ADD', width: 0.8, dash: 'dot'}},
      {type: 'line', x0: anio, x1: anio, y0: 0, y1: 1, yref: 'paper', line: {color: '#E24B4A', width: 1.5}}
    ]
  }, {displayModeBar: false, responsive: true});
}

// Mapa: estaciones de prueba coloreadas por la RECM de la KAN y estaciones de entrenamiento en gris
function dibujarMapa() {
  const est = D.por_anio[anio].estaciones;
  const capas = RELIEVE ? [{sourcetype: 'raster', source: ['https://tile.opentopomap.org/{z}/{x}/{y}.png'],
                            below: 'traces', opacity: 0.45, sourceattribution: 'OpenTopoMap (CC-BY-SA)'}] : [];
  Plotly.react('mapa', [
    {type: 'scattermap', lon: D.contexto.lon, lat: D.contexto.lat, mode: 'markers',
     marker: {size: 5, color: '#B4B2A9'}, hoverinfo: 'skip'},
    {type: 'scattermap', lon: est.map(e => e.lon), lat: est.map(e => e.lat), mode: 'markers',
     marker: {size: 13, color: est.map(e => e.recm), colorscale: 'Reds', cmin: 0, cmax: 25,
              colorbar: {title: {text: 'RECM KAN (mm)', side: 'right', font: {size: 11}}, thickness: 10, len: 0.5, x: 0.98}},
     customdata: est.map(e => [e.id, e.altitud, e.recm]),
     hovertemplate: 'Estación %{customdata[0]}<br>%{customdata[1]} m · RECM %{customdata[2]} mm<extra></extra>'}
  ], {
    map: {style: 'carto-positron', center: {lon: -76.55, lat: 3.95}, zoom: 7.1, layers: capas},
    margin: {l: 0, r: 0, t: 0, b: 0}, showlegend: false, font: fuente, paper_bgcolor: 'rgba(0,0,0,0)'
  }, {displayModeBar: false, responsive: true});
  const mapa = document.getElementById('mapa');
  mapa.removeAllListeners('plotly_hover');
  mapa.on('plotly_hover', ev => {
    const p = ev.points.find(q => q.curveNumber === 1);
    if (p) dibujarCurva(D.por_anio[anio].estaciones[p.pointIndex]);
  });
}

// Curva de la estación señalada: KAN, curva suavizada de la estación y satélite
function dibujarCurva(e) {
  actual = e;
  // El eje x son las etiquetas de cada péntada; en el eje se muestra el mes donde empieza
  const x = D.por_anio[anio].pentadas;
  const serie = (y, nombre, color, ancho, guion) => ({x: x, y: y, name: nombre, type: 'scatter', mode: 'lines',
    line: {color: color, width: ancho, dash: guion}, hovertemplate: nombre + ': %{y:.1f} mm<extra></extra>'});
  Plotly.react('curva', [
    serie(e.estacion, 'Estación', '#888780', 1.8, 'solid'),
    serie(e.satelite, 'Satélite', '#EF9F27', 1.4, 'dot'),
    serie(e.kan, 'KAN', '#E24B4A', 2.4, 'solid')
  ], {
    title: {text: 'Estación ' + e.id + ' · ' + e.altitud + ' m · ' + anio + ' · RECM ' + e.recm.toString().replace('.', ',') + ' mm',
            font: {size: 13}, x: 0.02},
    hovermode: 'x unified', margin: {l: 48, r: 12, t: 44, b: 44}, font: fuente,
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    xaxis: {type: 'category', showgrid: false, tickmode: 'array', tickvals: x.filter((_, i) => i % 6 === 0),
            ticktext: ['ene', 'feb', 'mar', 'abr', 'may', 'jun', 'jul', 'ago', 'sep', 'oct', 'nov', 'dic']},
    yaxis: {title: {text: 'Precipitación (mm por péntada)'}, rangemode: 'tozero', gridcolor: '#f1efe8'},
    legend: {orientation: 'h', y: -0.16, x: 0}
  }, {displayModeBar: false, responsive: true});
}

function elegirAnio(a) {
  anio = a;
  document.getElementById('anio').value = a;
  document.getElementById('etiqueta').textContent = a + ' · ' + faseAnio(a);
  dibujarOni();
  dibujarMapa();
  const est = D.por_anio[anio].estaciones;
  const misma = actual && est.find(e => e.id === actual.id);
  dibujarCurva(misma || est.reduce((m, e) => e.recm > m.recm ? e : m, est[0]));
}

document.getElementById('anio').addEventListener('input', ev => elegirAnio(+ev.target.value));
elegirAnio(1998);
document.getElementById('oni').on('plotly_click', ev => {
  const a = ev.points[0].x;
  if (a >= 1998 && a <= 2003) elegirAnio(a);
});
</script>
"""

pagina = (HTML.replace('PLOTLYJS', plotly.offline.get_plotlyjs())
              .replace('DATOS', json.dumps(datos))
              .replace('RELIEVE_ACTIVO', 'true' if relieve else 'false'))
components.html(pagina, height=740, scrolling=False)
st.caption('Curvas por péntada en mm. La estación es su curva suavizada; la KAN es la predicción del modelo elegido en la '
           'validación cruzada. Relieve: OpenTopoMap. Fase ENSO según el ONI promedio del año (umbral ±0,5).')
