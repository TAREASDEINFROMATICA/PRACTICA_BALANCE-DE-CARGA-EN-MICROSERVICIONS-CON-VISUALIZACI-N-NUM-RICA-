const EXACT = [120, 100, 80];
const CASES = {
  ideal: {
    name: 'Caso ideal',
    intro: 'Representa un sistema estable, con distribución de tráfico equilibrada, dependencias moderadas y una matriz bien condicionada. Es el escenario más favorable para observar una convergencia rápida y una interpretación clara del reparto de carga.',
    A: [[4,-1,-1],[-1,5,-1],[-1,-1,6]],
    b: [300,300,260]
  },
  stress: {
    name: 'Caso bajo estrés',
    intro: 'Representa un momento de alta demanda: muchos usuarios concurrentes, reintentos, cascadas de llamadas y mayor interacción entre microservicios. Los coeficientes crecen en magnitud y el problema se vuelve más exigente para los métodos iterativos.',
    A: [[100,-45,-30],[-45,90,-28],[-30,-28,85]],
    b: [5100,1360,400]
  },
  mal: {
    name: 'Caso mal condicionado',
    intro: 'Representa un entorno donde dos o más servicios se comportan casi igual ante el balanceador. Las ecuaciones resultan muy parecidas y los planos son casi paralelos. Esto vuelve al sistema numéricamente sensible y complica la convergencia.',
    A: [[1,0.99,0.98],[0.99,1,0.99],[0.98,0.99,1]],
    b: [297.4,298.0,296.6]
  }
};

function deepCopyMatrix(A){ return A.map(r=>r.slice()); }
function cloneVec(v){ return v.slice(); }
function fmt(x){
  if (typeof x !== 'number' || !isFinite(x)) return String(x);
  return Math.abs(x) >= 1000 || Math.abs(x) < 1e-3 && x !== 0 ? x.toExponential(4) : x.toFixed(6);
}
function round6(x){ return Math.round(x*1e6)/1e6; }
function vectorAdd(a,b){ return a.map((v,i)=>v+b[i]); }
function vectorSub(a,b){ return a.map((v,i)=>v-b[i]); }
function vectorScale(a,c){ return a.map(v=>v*c); }
function dot(a,b){ return a.reduce((s,v,i)=>s+v*b[i],0); }
function matVec(A,x){ return A.map(row => row.reduce((s,v,j)=>s+v*x[j],0)); }
function infNorm(v){ return Math.max(...v.map(n => Math.abs(n))); }
function euclid(v){ return Math.sqrt(dot(v,v)); }
function residualNorm(A,x,b){ return euclid(vectorSub(matVec(A,x),b)); }
function formatSystem(A,b){
  const lines = [];
  for(let i=0;i<3;i++){
    let s='';
    for(let j=0;j<3;j++){
      const a = A[i][j];
      const sign = a >= 0 ? (j===0 ? '' : ' + ') : ' - ';
      s += `${sign}${Math.abs(a)}x${j+1}`;
    }
    s += ` = ${b[i]}`;
    lines.push(s);
  }
  return lines;
}
function updateCaseUI(caseKey){
  const c = CASES[caseKey];
  document.querySelectorAll('.case-name').forEach(el=>el.textContent = c.name);
  document.querySelectorAll('.case-intro').forEach(el=>el.textContent = c.intro);
  const eqBox = document.getElementById('equations');
  if(eqBox){
    eqBox.innerHTML = formatSystem(c.A,c.b).map(line=>`<div class="eq-line">${line}</div>`).join('');
  }
  for(let i=0;i<3;i++){
    for(let j=0;j<3;j++){
      const inp = document.getElementById(`a${i}${j}`);
      if(inp) inp.value = c.A[i][j];
    }
    const binp = document.getElementById(`b${i}`);
    if(binp) binp.value = c.b[i];
  }
}
function readSystemFromInputs(){
  const A = [[],[],[]];
  const b = [];
  for(let i=0;i<3;i++){
    for(let j=0;j<3;j++) A[i][j] = parseFloat(document.getElementById(`a${i}${j}`).value);
    b[i] = parseFloat(document.getElementById(`b${i}`).value);
  }
  return {A,b};
}
function forwardSub(L,b, unit=false){
  const n = b.length, y = Array(n).fill(0);
  for(let i=0;i<n;i++){
    let s = b[i];
    for(let j=0;j<i;j++) s -= L[i][j]*y[j];
    y[i] = unit ? s : s/L[i][i];
  }
  return y;
}
function backwardSub(U,y, unit=false){
  const n = y.length, x = Array(n).fill(0);
  for(let i=n-1;i>=0;i--){
    let s = y[i];
    for(let j=i+1;j<n;j++) s -= U[i][j]*x[j];
    x[i] = unit ? s : s/U[i][i];
  }
  return x;
}
function luDoolittle(A){
  const n = A.length;
  const L = Array.from({length:n},()=>Array(n).fill(0));
  const U = Array.from({length:n},()=>Array(n).fill(0));
  for(let i=0;i<n;i++){
    for(let k=i;k<n;k++){
      let sum=0; for(let j=0;j<i;j++) sum += L[i][j]*U[j][k];
      U[i][k] = A[i][k]-sum;
    }
    L[i][i]=1;
    for(let k=i+1;k<n;k++){
      let sum=0; for(let j=0;j<i;j++) sum += L[k][j]*U[j][i];
      L[k][i] = (A[k][i]-sum)/U[i][i];
    }
  }
  return {L,U};
}
function luCrout(A){
  const n=A.length;
  const L=Array.from({length:n},()=>Array(n).fill(0));
  const U=Array.from({length:n},()=>Array(n).fill(0));
  for(let j=0;j<n;j++) U[j][j]=1;
  for(let j=0;j<n;j++){
    for(let i=j;i<n;i++){
      let sum=0; for(let k=0;k<j;k++) sum += L[i][k]*U[k][j];
      L[i][j] = A[i][j]-sum;
    }
    for(let i=j+1;i<n;i++){
      let sum=0; for(let k=0;k<j;k++) sum += L[j][k]*U[k][i];
      U[j][i] = (A[j][i]-sum)/L[j][j];
    }
  }
  return {L,U};
}
function cholesky(A){
  const n=A.length;
  const L=Array.from({length:n},()=>Array(n).fill(0));
  for(let i=0;i<n;i++){
    for(let j=0;j<=i;j++){
      let sum=0; for(let k=0;k<j;k++) sum += L[i][k]*L[j][k];
      if(i===j){
        const val = A[i][i]-sum;
        if(val <= 0) throw new Error('La matriz no es definida positiva, Cholesky no puede aplicarse directamente.');
        L[i][j] = Math.sqrt(val);
      } else {
        L[i][j] = (A[i][j]-sum)/L[j][j];
      }
    }
  }
  const U = transpose(L);
  return {L,U};
}
function transpose(A){ return A[0].map((_,j)=>A.map(r=>r[j])); }
function solveWithLU(A,b, variant){
  let fac;
  if(variant==='doolittle') fac = luDoolittle(A);
  else if(variant==='crout') fac = luCrout(A);
  else fac = cholesky(A);
  const unitL = variant==='doolittle';
  const unitU = variant==='crout';
  const y = forwardSub(fac.L,b, unitL);
  const x = backwardSub(fac.U,y, unitU);
  return { ...fac, y, x, residual: residualNorm(A,x,b) };
}
function jacobi(A,b,tol=1e-6,maxIter=200,x0=[0,0,0]){
  const n=b.length; let x=x0.slice(); let hist=[];
  for(let k=1;k<=maxIter;k++){
    const xnew=Array(n).fill(0);
    for(let i=0;i<n;i++){
      let s=b[i];
      for(let j=0;j<n;j++) if(j!==i) s -= A[i][j]*x[j];
      xnew[i] = s/A[i][i];
    }
    const err=infNorm(vectorSub(xnew,x));
    const exactErr = infNorm(vectorSub(xnew,EXACT));
    hist.push({k,x:xnew.slice(),err,exactErr,res:residualNorm(A,xnew,b)});
    x=xnew;
    if(err<tol) return {x, iterations:k, converged:true, history:hist};
  }
  return {x, iterations:maxIter, converged:false, history:hist};
}
function gaussSeidel(A,b,tol=1e-6,maxIter=200,x0=[0,0,0]){
  const n=b.length; let x=x0.slice(); let hist=[];
  for(let k=1;k<=maxIter;k++){
    const old=x.slice();
    for(let i=0;i<n;i++){
      let s=b[i];
      for(let j=0;j<n;j++) if(j!==i) s -= A[i][j]*x[j];
      x[i]=s/A[i][i];
    }
    const err=infNorm(vectorSub(x,old));
    hist.push({k,x:x.slice(),err,exactErr:infNorm(vectorSub(x,EXACT)),res:residualNorm(A,x,b)});
    if(err<tol) return {x, iterations:k, converged:true, history:hist};
  }
  return {x, iterations:maxIter, converged:false, history:hist};
}
function sor(A,b,w=1.15,tol=1e-6,maxIter=200,x0=[0,0,0]){
  const n=b.length; let x=x0.slice(); let hist=[];
  for(let k=1;k<=maxIter;k++){
    const old=x.slice();
    for(let i=0;i<n;i++){
      let s1=0,s2=0;
      for(let j=0;j<i;j++) s1 += A[i][j]*x[j];
      for(let j=i+1;j<n;j++) s2 += A[i][j]*old[j];
      const gs=(b[i]-s1-s2)/A[i][i];
      x[i]=(1-w)*old[i]+w*gs;
    }
    const err=infNorm(vectorSub(x,old));
    hist.push({k,x:x.slice(),err,exactErr:infNorm(vectorSub(x,EXACT)),res:residualNorm(A,x,b)});
    if(err<tol) return {x, iterations:k, converged:true, history:hist};
  }
  return {x, iterations:maxIter, converged:false, history:hist};
}
function pcg(A,b,tol=1e-6,maxIter=100,x0=[0,0,0]){
  const n=b.length;
  let x=x0.slice();
  let r=vectorSub(b,matVec(A,x));
  const Mdiag=A.map((row,i)=>row[i]);
  let z=r.map((ri,i)=>ri/Mdiag[i]);
  let p=z.slice();
  let rzOld=dot(r,z);
  let hist=[];
  for(let k=1;k<=maxIter;k++){
    const Ap=matVec(A,p);
    const alpha=rzOld/dot(p,Ap);
    x=vectorAdd(x, vectorScale(p,alpha));
    r=vectorSub(r, vectorScale(Ap,alpha));
    const res=euclid(r);
    hist.push({k,x:x.slice(),err:res,exactErr:infNorm(vectorSub(x,EXACT)),res});
    if(res<tol) return {x, iterations:k, converged:true, history:hist};
    z=r.map((ri,i)=>ri/Mdiag[i]);
    const rzNew=dot(r,z);
    const beta=rzNew/rzOld;
    p=vectorAdd(z, vectorScale(p,beta));
    rzOld=rzNew;
  }
  return {x, iterations:maxIter, converged:false, history:hist};
}

function matrixHtml(M){
  return `<pre>${M.map(r=>r.map(v=>fmt(v)).join('   ')).join('\n')}</pre>`;
}
function vectorHtml(v){ return `<pre>${v.map(fmt).join('\n')}</pre>`; }
function buildIterRows(hist, n=6){
  const rows = hist.slice(0,n).map(h=>`<tr><td>${h.k}</td><td>${fmt(h.x[0])}</td><td>${fmt(h.x[1])}</td><td>${fmt(h.x[2])}</td><td>${fmt(h.err)}</td><td>${fmt(h.res)}</td></tr>`).join('');
  return rows || '<tr><td colspan="6">No aplica</td></tr>';
}

function drawLineChart(canvas, datasets, labels, title){
  const ctx = canvas.getContext('2d');
  const W = canvas.width = canvas.clientWidth * (window.devicePixelRatio || 1);
  const H = canvas.height = canvas.clientHeight * (window.devicePixelRatio || 1);
  const ratio = (window.devicePixelRatio || 1);
  ctx.setTransform(ratio,0,0,ratio,0,0);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  ctx.clearRect(0,0,w,h);
  const pad = {l:50,r:20,t:35,b:45};
  const plotW = w-pad.l-pad.r, plotH = h-pad.t-pad.b;
  const allValues = datasets.flatMap(d => d.data).filter(v => isFinite(v));
  const maxY = Math.max(...allValues, 1);
  const minY = 0;
  ctx.strokeStyle='#cbd5e1'; ctx.lineWidth=1;
  for(let i=0;i<=5;i++){
    const y = pad.t + plotH*i/5;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(w-pad.r,y); ctx.stroke();
    const val = maxY - (maxY-minY)*i/5;
    ctx.fillStyle='#475569'; ctx.font='12px Arial'; ctx.fillText(val.toFixed(2),5,y+4);
  }
  ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,h-pad.b); ctx.lineTo(w-pad.r,h-pad.b); ctx.strokeStyle='#334155'; ctx.lineWidth=1.5; ctx.stroke();
  const colors=['#2563eb','#dc2626','#059669','#d97706','#7c3aed'];
  datasets.forEach((ds,idx)=>{
    ctx.strokeStyle=colors[idx%colors.length]; ctx.lineWidth=2; ctx.beginPath();
    ds.data.forEach((v,i)=>{
      const x = pad.l + (labels.length===1?0:i*(plotW/(labels.length-1)));
      const y = pad.t + plotH*(1-(v-minY)/(maxY-minY || 1));
      if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    });
    ctx.stroke();
    ds.data.forEach((v,i)=>{
      const x = pad.l + (labels.length===1?0:i*(plotW/(labels.length-1)));
      const y = pad.t + plotH*(1-(v-minY)/(maxY-minY || 1));
      ctx.fillStyle=colors[idx%colors.length]; ctx.beginPath(); ctx.arc(x,y,3,0,Math.PI*2); ctx.fill();
    });
  });
  ctx.fillStyle='#0f172a'; ctx.font='bold 15px Arial'; ctx.fillText(title, pad.l, 20);
  ctx.font='12px Arial';
  labels.forEach((lab,i)=>{
    if(labels.length>25 && i%Math.ceil(labels.length/10)!==0 && i!==labels.length-1) return;
    const x = pad.l + (labels.length===1?0:i*(plotW/(labels.length-1)));
    ctx.fillStyle='#475569'; ctx.fillText(String(lab), x-5, h-18);
  });
  datasets.forEach((ds,idx)=>{ ctx.fillStyle=colors[idx%colors.length]; ctx.fillRect(w-160, 16+idx*18, 12, 12); ctx.fillStyle='#334155'; ctx.fillText(ds.label, w-142, 26+idx*18); });
}
function drawBarChart(canvas, datasets, categories, title){
  const ctx = canvas.getContext('2d');
  const W = canvas.width = canvas.clientWidth * (window.devicePixelRatio || 1);
  const H = canvas.height = canvas.clientHeight * (window.devicePixelRatio || 1);
  const ratio = (window.devicePixelRatio || 1);
  ctx.setTransform(ratio,0,0,ratio,0,0);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  ctx.clearRect(0,0,w,h);
  const pad={l:50,r:20,t:35,b:45};
  const plotW=w-pad.l-pad.r, plotH=h-pad.t-pad.b;
  const maxY=Math.max(...datasets.flatMap(d=>d.data), 1);
  ctx.strokeStyle='#334155'; ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,h-pad.b); ctx.lineTo(w-pad.r,h-pad.b); ctx.stroke();
  ctx.strokeStyle='#e2e8f0';
  for(let i=0;i<=5;i++){
    const y = pad.t + plotH*i/5; ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(w-pad.r,y); ctx.stroke();
    const val = maxY - maxY*i/5; ctx.fillStyle='#475569'; ctx.font='12px Arial'; ctx.fillText(val.toFixed(2),5,y+4);
  }
  const colors=['#2563eb','#dc2626','#059669','#d97706','#7c3aed'];
  const groups=categories.length, dsCount=datasets.length;
  const groupW=plotW/groups;
  const barW=Math.min(28,(groupW*0.75)/dsCount);
  categories.forEach((cat,i)=>{
    const start = pad.l + i*groupW + groupW*0.12;
    datasets.forEach((ds,j)=>{
      const x = start + j*barW;
      const bh = (ds.data[i]/maxY)*plotH;
      const y = h-pad.b-bh;
      ctx.fillStyle = colors[j%colors.length];
      ctx.fillRect(x,y,barW-4,bh);
    });
    ctx.fillStyle='#475569'; ctx.fillText(cat, start + (dsCount*barW)/2 - 10, h-18);
  });
  ctx.fillStyle='#0f172a'; ctx.font='bold 15px Arial'; ctx.fillText(title, pad.l, 20);
  datasets.forEach((ds,idx)=>{ ctx.fillStyle=colors[idx%colors.length]; ctx.fillRect(w-160, 16+idx*18, 12, 12); ctx.fillStyle='#334155'; ctx.fillText(ds.label, w-142, 26+idx*18); });
}

function renderLU() {
  const { A, b } = readSystemFromInputs();

  const variants = ["doolittle", "crout", "cholesky"];
  const names = {
    doolittle: "Doolittle",
    crout: "Crout",
    cholesky: "Cholesky"
  };

  let html = "";
  const residuals = [];
  const solutionSets = [];

  variants.forEach(v => {
    try {
      const res = solveWithLU(deepCopyMatrix(A), cloneVec(b), v);

      residuals.push(res.residual);
      solutionSets.push(res.x);

      html += `
        <div class="card">
          <h3>${names[v]}</h3>

          <div class="grid-2">
            <div>
              <strong>Matriz L</strong>
              ${matrixHtml(res.L)}
            </div>
            <div>
              <strong>Matriz U</strong>
              ${matrixHtml(res.U)}
            </div>
          </div>

          <p><strong>Vector intermedio y:</strong></p>
          ${vectorHtml(res.y)}

          <p>
            <strong>Solución final:</strong><br>
            x1 = ${fmt(res.x[0])}<br>
            x2 = ${fmt(res.x[1])}<br>
            x3 = ${fmt(res.x[2])}
          </p>

          <p><strong>Norma del residuo:</strong> ${fmt(res.residual)}</p>

          <div class="result-box" style="margin-top:12px">
            <strong>Análisis:</strong>
            <p>${analysisLU(v, res.residual)}</p>
          </div>

          <div class="result-box" style="margin-top:12px">
            <strong>Conclusión:</strong>
            <p>${conclusionLU(v, res.residual)}</p>
          </div>
        </div>
      `;
    } catch (err) {
      residuals.push(NaN);
      solutionSets.push([NaN, NaN, NaN]);

      html += `
        <div class="card">
          <h3>${names[v]}</h3>
          <p>${err.message}</p>
        </div>
      `;
    }
  });

  document.getElementById("results").innerHTML = html;

  document.getElementById("final-solution").innerHTML = `
    <div class="grid-3" style="margin-top:12px">
      <div class="kpi">
        <span class="small">Solución exacta esperada x1</span>
        <strong>${fmt(EXACT[0])}</strong>
      </div>
      <div class="kpi">
        <span class="small">Solución exacta esperada x2</span>
        <strong>${fmt(EXACT[1])}</strong>
      </div>
      <div class="kpi">
        <span class="small">Solución exacta esperada x3</span>
        <strong>${fmt(EXACT[2])}</strong>
      </div>
    </div>
  `;

  document.getElementById("analysis-text").innerHTML =
    "La factorización LU permite resolver el sistema de forma directa, sin depender de iteraciones sucesivas. " +
    "Esto la convierte en un método de referencia muy útil dentro del problema de balance de carga en microservicios, " +
    "porque proporciona una solución estable con la cual se pueden comparar Jacobi, Gauss-Seidel, SOR y Gradiente Conjugado Precondicionado. " +
    "Además, las variantes Doolittle, Crout y Cholesky muestran distintas formas de descomponer la matriz, pero todas buscan llegar a la misma solución del sistema.";

  document.getElementById("conclusion-text").innerHTML =
    "En conclusión, el método LU es apropiado cuando se desea una solución directa, precisa y verificable. " +
    "Dentro de este proyecto, LU sirve como base de contraste frente a los métodos iterativos. " +
    "Su principal ventaja es que no necesita iteraciones para converger, por lo que reduce el riesgo de inestabilidad numérica en escenarios favorables y permite comprobar si las aproximaciones obtenidas por otros métodos son correctas.";

  document.getElementById("conclusion-extra").innerHTML =
    "Doolittle y Crout son dos formas equivalentes de factorizar la matriz, diferenciándose en la ubicación de la diagonal unitaria. " +
    "Cholesky, en cambio, resulta más eficiente cuando la matriz cumple la condición de ser simétrica definida positiva. " +
    "Por eso, en un análisis comparativo, Doolittle y Crout son útiles como opciones generales, mientras que Cholesky destaca cuando la estructura de la matriz lo permite. " +
    "Esto demuestra que no existe una única variante mejor para todos los casos, sino que la elección depende de las propiedades del sistema.";

  document.getElementById("conclusion-microservices").innerHTML =
    "Aplicado al balance de carga en microservicios, LU permite interpretar con claridad cómo se reparte la carga entre autenticación, catálogo y pagos. " +
    "Si la solución obtenida coincide con el comportamiento esperado, se puede afirmar que el sistema modela correctamente la distribución de solicitudes. " +
    "Además, en escenarios como el caso ideal, bajo estrés o mal condicionado, LU ayuda a verificar si las diferencias observadas en los métodos iterativos se deben al algoritmo o al condicionamiento del sistema. " +
    "Por eso, LU no solo resuelve el problema, sino que también cumple una función de validación y referencia en el análisis del rendimiento del sistema.";

  drawBarChart(
    document.getElementById("chart1"),
    [
      { label: "Exacta", data: EXACT },
      { label: "Doolittle", data: solutionSets[0] },
      { label: "Crout", data: solutionSets[1] },
      { label: "Cholesky", data: solutionSets[2] }
    ],
    ["x1", "x2", "x3"],
    "Comparación de soluciones obtenidas"
  );

  drawBarChart(
    document.getElementById("chart2"),
    [
      {
        label: "Residuo",
        data: residuals.map(v => isFinite(v) ? v : 0)
      }
    ],
    ["Doolittle", "Crout", "Cholesky"],
    "Norma del residuo por variante LU"
  );
}
function analysisLU(v,res){
  const base={doolittle:'Doolittle separa la matriz en una triangular inferior con diagonal unitaria y una triangular superior completa.',crout:'Crout usa el esquema inverso: la diagonal unitaria está en U y la información principal queda en L.',cholesky:'Cholesky aprovecha que la matriz es simétrica definida positiva para factorizar A = L·Lᵀ de forma más eficiente.'};
  return `${base[v]} En este escenario, el residuo obtenido es ${fmt(res)}, lo que muestra que la solución directa reproduce el balance esperado con gran precisión.`;
}
function conclusionLU(v,res){
  return `La variante ${v==='doolittle'?'Doolittle':v==='crout'?'Crout':'Cholesky'} es adecuada como referencia exacta. Al no depender de iteraciones, sirve para validar los métodos iterativos y comparar estabilidad en los tres casos.`;
}
function renderIterative(method){
  const {A,b}=readSystemFromInputs();
  const tol=1e-6; const maxIter=method==='jacobi'?120:(method==='pcg'?30:120);
  let res;
  if(method==='jacobi') res = jacobi(A,b,tol,maxIter);
  if(method==='gs') res = gaussSeidel(A,b,tol,maxIter);
  if(method==='sor') {
    const omega = parseFloat(document.getElementById('omega').value);
    res = sor(A,b,omega,tol,maxIter);
  }
  if(method==='pcg') res = pcg(A,b,tol,maxIter);

  document.getElementById('final-solution').innerHTML = `
    <div class="grid-3">
      <div class="kpi"><span class="small">x1</span><strong>${fmt(res.x[0])}</strong></div>
      <div class="kpi"><span class="small">x2</span><strong>${fmt(res.x[1])}</strong></div>
      <div class="kpi"><span class="small">x3</span><strong>${fmt(res.x[2])}</strong></div>
    </div>
    <div class="grid-3" style="margin-top:12px">
      <div class="kpi"><span class="small">Iteraciones</span><strong>${res.iterations}</strong></div>
      <div class="kpi"><span class="small">Convergencia</span><strong>${res.converged?'Sí':'No'}</strong></div>
      <div class="kpi"><span class="small">Residuo final</span><strong>${fmt(residualNorm(A,res.x,b))}</strong></div>
    </div>`;
  document.getElementById('iter-body').innerHTML = buildIterRows(res.history,8);
  document.getElementById('analysis-text').innerHTML = analysisIterative(method,res);
  document.getElementById('conclusion-text').innerHTML = conclusionIterative(method,res);
  drawLineChart(document.getElementById('chart1'), [
    {label:'Error respecto a solución exacta', data: res.history.map(h=>h.exactErr)},
    {label:'Residuo', data: res.history.map(h=>h.res)}
  ], res.history.map(h=>h.k), 'Convergencia por iteración');
  drawBarChart(document.getElementById('chart2'), [
    {label:'Exacta', data: EXACT},
    {label:'Aproximada', data: res.x}
  ], ['x1','x2','x3'], 'Solución final vs solución exacta');
}
function analysisIterative(method,res){
  const names={jacobi:'Jacobi', gs:'Gauss-Seidel', sor:'SOR', pcg:'Gradiente Conjugado Precondicionado'};
  const base={
    jacobi:'Jacobi actualiza todas las variables utilizando únicamente los valores de la iteración previa, por eso suele ser el método más conservador y, a la vez, el más lento.',
    gs:'Gauss-Seidel usa de inmediato los nuevos valores calculados dentro de la misma iteración. Eso normalmente acelera la convergencia frente a Jacobi.',
    sor:'SOR añade un parámetro de relajación ω para adelantar la aproximación y reducir el número de iteraciones cuando ω se elige adecuadamente.',
    pcg:'El Gradiente Conjugado Precondicionado trabaja con direcciones conjugadas y usa un precondicionador diagonal. En matrices simétricas definidas positivas suele ser muy eficiente.'
  };
  return `${base[method]} En la ejecución actual, el método ${res.converged?'sí logró converger':'no alcanzó la tolerancia establecida'} y necesitó <strong>${res.iterations}</strong> iteraciones. La solución aproximada obtenida fue (${fmt(res.x[0])}, ${fmt(res.x[1])}, ${fmt(res.x[2])}) con residuo final ${fmt(res.history.length ? res.history[res.history.length-1].res : residualNorm(readSystemFromInputs().A,res.x,readSystemFromInputs().b))}. Esto permite evaluar si el escenario escogido es amigable, exigente o numéricamente problemático.`;
}
function conclusionIterative(method,res){
  const m = {jacobi:'Jacobi', gs:'Gauss-Seidel', sor:'SOR', pcg:'Gradiente Conjugado Precondicionado'}[method];
  return `Conclusión: ${m} ${res.converged?`es viable para este escenario, porque alcanzó la tolerancia y produjo una solución cercana a la exacta en ${res.iterations} iteraciones.`:`muestra dificultades para este escenario dentro del máximo de iteraciones planteado, lo que evidencia sensibilidad al condicionamiento o una velocidad de convergencia insuficiente.`} Esto lo convierte en una pieza útil para comparar estabilidad y costo computacional frente a los demás métodos.`;
}
function renderPage(){
  const page = document.body.dataset.page;
  const caseSelect = document.getElementById('caseSelect');
  if(caseSelect){
    updateCaseUI(caseSelect.value);
    caseSelect.addEventListener('change', ()=>updateCaseUI(caseSelect.value));
  }
  const fillBtn = document.getElementById('fillDefault');
  if(fillBtn) fillBtn.addEventListener('click', ()=>updateCaseUI(caseSelect.value));
  if(page==='lu'){
    document.getElementById('solveBtn').addEventListener('click', renderLU);
    renderLU();
  }
  if(page==='jacobi'){
    document.getElementById('solveBtn').addEventListener('click', ()=>renderIterative('jacobi'));
    renderIterative('jacobi');
  }
  if(page==='gauss'){
    document.getElementById('solveBtn').addEventListener('click', ()=>renderIterative('gs'));
    renderIterative('gs');
  }
  if(page==='sor'){
    document.getElementById('solveBtn').addEventListener('click', ()=>renderIterative('sor'));
    renderIterative('sor');
  }
  if(page==='pcg'){
    document.getElementById('solveBtn').addEventListener('click', ()=>renderIterative('pcg'));
    renderIterative('pcg');
  }
  window.addEventListener('resize', ()=>{
    if(page==='lu') renderLU();
    if(page==='jacobi') renderIterative('jacobi');
    if(page==='gauss') renderIterative('gs');
    if(page==='sor') renderIterative('sor');
    if(page==='pcg') renderIterative('pcg');
  });
}
window.addEventListener('DOMContentLoaded', renderPage);
