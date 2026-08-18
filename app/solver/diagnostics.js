/** Diagnostics shared by the sixfold overview and future field solvers. */

export function sampleRadialField(fieldFn, { radius, maxRadius, arms = 6, radialSamples = 42 }) {
  const result = [];
  for (let i = 0; i < arms; i++) {
    const theta = i * 2 * Math.PI / arms;
    const samples = [];
    for (let j = 0; j < radialSamples; j++) {
      const r = radius + (maxRadius - radius) * j / (radialSamples - 1);
      samples.push({ r, value: fieldFn(r, theta) });
    }
    result.push({ theta, samples });
  }
  return result;
}

/**
 * Dense angular Fourier decomposition. Do not use only six crystallographic rays:
 * m=6 aliases m=0 on those samples. Dense ring sampling avoids that degeneracy.
 */
export function angularModeSpectrum(fieldFn, { r, angularSamples = 72, modes = [0, 1, 2, 3, 6, 12] }) {
  const samples = [];
  for (let k = 0; k < angularSamples; k++) {
    const theta = 2 * Math.PI * k / angularSamples;
    samples.push({ theta, value: fieldFn(r, theta) });
  }
  const mean = samples.reduce((s, p) => s + p.value, 0) / samples.length;
  return modes.map(m => {
    if (m === 0) return { m, a: mean, b: 0, amplitude: Math.abs(mean), phase: 0 };
    let a = 0, b = 0;
    for (const p of samples) {
      const residual = p.value - mean;
      a += residual * Math.cos(m * p.theta);
      b += residual * Math.sin(m * p.theta);
    }
    a *= 2 / samples.length;
    b *= 2 / samples.length;
    return { m, a, b, amplitude: Math.hypot(a, b), phase: Math.atan2(-b, a) / m };
  });
}

export function sixfoldDeviation(values) {
  const mean = values.reduce((a, b) => a + b, 0) / Math.max(1, values.length);
  return { mean, residuals: values.map(v => v - mean) };
}
