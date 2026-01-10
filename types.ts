
export interface BattenInputs {
  testWeight: number;
  testLength: number;
  self14: number;
  self12: number;
  self34: number;
  weighted14: number;
  weighted12: number;
  weighted34: number;
}

export interface BattenResults {
  frontPercent: number;
  backPercent: number;
  camberPercent: number;
  averageEi: number;
  deflection: number;
  draftPosition: number;
}

export interface SavedProfile {
  id: string;
  name: string;
  inputs: BattenInputs;
  date: string;
}

export interface BattenSegment {
  length: number;
  ei: number;
}

export interface FEMCalibrationInputs {
  weightKg: number;
  lengthMm: number;
  deflectionsMm: number[];
}

export interface FEMCalibrationResults {
  eiSegments: number[];
  eiAverage: number;
  error: number;
}
