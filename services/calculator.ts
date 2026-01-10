
import { BattenInputs, BattenResults, BattenSegment } from '../types';

/**
 * Calcola l'EI equivalente di una stecca composta da N segmenti di cui si conosce EI e L.
 * Basato sull'integrale della curvatura per deflessione al centro.
 */
export const calculateEquivalentEI = (segments: BattenSegment[]): { totalLength: number, eqEi: number } => {
  const totalLength = segments.reduce((acc, s) => acc + s.length, 0);
  if (totalLength <= 0) return { totalLength: 0, eqEi: 0 };

  const halfL = totalLength / 2;
  let integral = 0;
  let currentPos = 0;

  // Integriamo x^2 / EI(x) da 0 a L/2
  for (const seg of segments) {
    const start = currentPos;
    const end = Math.min(currentPos + seg.length, halfL);
    
    if (end > start && seg.ei > 0) {
      const ei_mm = seg.ei * 1e6; // Converti Nm2 in Nmm2
      integral += (Math.pow(end, 3) - Math.pow(start, 3)) / (3 * ei_mm);
    }
    currentPos += seg.length;
    if (currentPos >= halfL) break;
  }

  // Formula EI_eq = L^3 / (48 * delta) dove delta = F * integral / (??) 
  // In realtà EI_eq = L^3 / (24 * integral_0_L/2) per trave simmetrica
  const eqEi = (Math.pow(totalLength, 3) / (24 * integral)) / 1e6;

  return { 
    totalLength, 
    eqEi: isFinite(eqEi) ? Number(eqEi.toFixed(3)) : 0 
  };
};

export const calculateBattenBehavior = (inputs: BattenInputs): BattenResults => {
  const { testWeight, testLength, self14, self12, self34, weighted14, weighted12, weighted34 } = inputs;
  const d14 = Math.max(0, weighted14 - self14);
  const d12 = Math.max(0.1, weighted12 - self12);
  const d34 = Math.max(0, weighted34 - self34);
  const Lm = testLength / 1000;
  const F_newtons = testWeight * 9.80665;
  const camberPercent = testLength > 0 ? (d12 / testLength) * 100 : 0;
  const ei = d12 > 0 ? (F_newtons * Math.pow(Lm, 3)) / (48 * (d12 / 1000)) : 0;
  const frontPercent = (d14 / d12) * 100;
  const backPercent = (d34 / d12) * 100;
  const draftPos = 50 - ((d14 - d34) / d12 * 10);

  return {
    frontPercent: Number(frontPercent.toFixed(1)),
    backPercent: Number(backPercent.toFixed(1)),
    camberPercent: Number(camberPercent.toFixed(2)),
    averageEi: Number(ei.toFixed(3)),
    deflection: d12,
    draftPosition: Number(draftPos.toFixed(1))
  };
};

export const generateCurvePoints = (inputs: BattenInputs) => {
  const { testLength, self14, self12, self34, weighted14, weighted12, weighted34 } = inputs;
  const d14 = Math.max(0, weighted14 - self14);
  const d12 = Math.max(0, weighted12 - self12);
  const d34 = Math.max(0, weighted34 - self34);
  const points = [];
  const steps = 60;
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const x = t * testLength;
    const l1 = (t - 0) * (t - 0.5) * (t - 0.75) * (t - 1.0) / ( (0.25 - 0) * (0.25 - 0.5) * (0.25 - 0.75) * (0.25 - 1.0) );
    const l2 = (t - 0) * (t - 0.25) * (t - 0.75) * (t - 1.0) / ( (0.5 - 0) * (0.5 - 0.25) * (0.5 - 0.75) * (0.5 - 1.0) );
    const l3 = (t - 0) * (t - 0.25) * (t - 0.5) * (t - 1.0) / ( (0.75 - 0) * (0.75 - 0.25) * (0.75 - 0.5) * (0.75 - 1.0) );
    const y = (d14 * l1) + (d12 * l2) + (d34 * l3);
    points.push({ x: Math.round(x), y: Number(y.toFixed(2)) });
  }
  return points;
};
