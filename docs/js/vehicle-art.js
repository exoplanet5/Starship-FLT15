// Cartoon of the launch vehicle, drawn nozzle-down on a square canvas from meta.stack: segment lengths and widths in metres, with the width doubled so the
// shape reads at sprite size. The stack shortens as stages separate, so each state has its own segments and scale. A state may carry side
// boosters, {count, len, w, burn}: two are drawn in profile beside the core, with their own flames while burn is set.
const COL = {outline: '#0a0d12', body: ['#8d97a2', '#f6f8fa', '#d7dde3', '#7c8792'], dark: ['#27303a', '#59656f', '#39434d', '#1a2027'], band: '#c8372d', gold: ['#8a6a1f', '#ecc75a', '#c79c36', '#73581a']};
const WIDTH_BOOST = 2.0;   // default; a stack may set its own width_boost, 1 for true proportions

// Missions without meta.stack get a plain rocket: the full stack under power, then the payload.
export const genericStack = (meta) => ({
  total_m: 30, burn_end: meta.insertion_met,
  segments: {core: {len: 18, w: 3, kind: 'solid', band: true}, upper: {len: 4, w: 3, kind: 'liquid'}, fairing: {len: 8, w: 3.6, kind: 'fairing'}, payload: {len: 3, w: 2, kind: 'payload'}},
  states: [{until: meta.insertion_met, segs: ['core', 'upper', 'fairing'], zoom: 1, burn: 'liquid'}, {segs: ['payload'], zoom: 3.2, burn: null}],
});

export const vehicleStateAt = (met, stack, poweredAt) => {
  const found = stack.states.findIndex((s) => s.until !== undefined && met < s.until);
  return {state: found < 0 ? stack.states.length - 1 : found, burning: poweredAt(met) && met <= stack.burn_end};
};

function shade(g, x, w, stops) {
  const grad = g.createLinearGradient(x - w / 2, 0, x + w / 2, 0);
  grad.addColorStop(0, stops[0]); grad.addColorStop(.32, stops[1]); grad.addColorStop(.68, stops[2]); grad.addColorStop(1, stops[3]);
  return grad;
}
function flame(g, x, y, w, len, kind) {
  const solid = kind === 'solid';
  const plume = (wid, l, stops) => {
    const grad = g.createLinearGradient(0, y, 0, y + l);
    stops.forEach(([p, c]) => grad.addColorStop(p, c));
    g.fillStyle = grad; g.beginPath(); g.moveTo(x - wid / 2, y);
    g.bezierCurveTo(x - wid * .62, y + l * .3, x - wid * .22, y + l * .72, x, y + l);
    g.bezierCurveTo(x + wid * .22, y + l * .72, x + wid * .62, y + l * .3, x + wid / 2, y); g.closePath(); g.fill();
  };
  if (solid) {
    plume(w, len, [[0, 'rgba(255,190,80,.95)'], [.45, 'rgba(255,120,30,.72)'], [1, 'rgba(255,70,10,0)']]);
    plume(w * .5, len * .62, [[0, 'rgba(255,255,240,1)'], [.6, 'rgba(255,225,130,.85)'], [1, 'rgba(255,190,80,0)']]);
  } else {
    plume(w, len, [[0, 'rgba(150,205,255,.85)'], [1, 'rgba(90,160,255,0)']]);
    plume(w * .45, len * .6, [[0, 'rgba(240,248,255,1)'], [1, 'rgba(170,215,255,0)']]);
  }
}

// Starship and Super Heavy, kinds 'starship' and 'superheavy', in profile: stainless hulls with only the parts that identify them. Each is drawn from its
// base (x, bot) upward over length l at width w; heights are fractions of l and spans fractions of w, so the shapes keep their proportions at any width_boost.
const STEEL = ['#858f9a', '#f4f6f8', '#cdd3d9', '#747e89'], TILE = '#232830', DARK = '#2b3139';
const poly = (g, pts) => { g.beginPath(); pts.forEach(([px, py], i) => (i ? g.lineTo(px, py) : g.moveTo(px, py))); g.closePath(); g.fill(); g.stroke(); };
const OWN_BASE = new Set(['starship', 'superheavy']);   // these draw their own engine section, so the stack adds no nozzle

// Ship: barrel and blunt ogive nose (29 % of the length), tiled windward side on the left, large aft flaps and small forward flaps set on the nose.
function starship(g, x, bot, w, l) {
  const top = bot - l, nose = l * .29, sh = top + nose, L = x - w / 2, R = x + w / 2;
  for (const s of [-1, 1]) {
    const e = x + s * w / 2;
    g.fillStyle = s < 0 ? TILE : STEEL[3];
    poly(g, [[e, bot - l * .29], [e + s * w * .48, bot - l * .16], [e + s * w * .48, bot - l * .015], [e, bot - l * .015]]);                 // aft flap
    poly(g, [[x + s * w * .15, top + l * .07], [e + s * w * .2, top + l * .185], [e + s * w * .2, top + l * .25], [x + s * w * .3, top + l * .275]]);   // forward flap
  }
  const body = new Path2D();
  body.moveTo(L, bot); body.lineTo(L, sh); body.bezierCurveTo(L, sh - nose * .58, x - w * .17, top, x, top);
  body.bezierCurveTo(x + w * .17, top, R, sh - nose * .58, R, sh); body.lineTo(R, bot); body.closePath();
  g.fillStyle = shade(g, x, w, STEEL); g.fill(body);
  g.save(); g.clip(body); g.fillStyle = TILE; g.fillRect(L, top, w * .46, l); g.restore();
  g.stroke(body);
}

// Booster: plain barrel with the hot-stage ring on top, the engine section at the base, and grid fins edge-on below the ring. The V3 booster has three
// fins; in this profile the third is behind the barrel.
function superHeavy(g, x, bot, w, l) {
  const top = bot - l, L = x - w / 2, fin = top + l * .1, th = l * .026;
  g.fillStyle = DARK;
  for (const s of [-1, 1]) g.fillRect(x + s * w / 2, fin, s * w * .42, th), g.strokeRect(x + s * w / 2, fin, s * w * .42, th);
  g.fillStyle = shade(g, x, w, STEEL); g.fillRect(L, top, w, l);
  g.fillStyle = DARK; g.fillRect(L, top, w, l * .04); g.fillRect(L, bot - l * .03, w, l * .03);
  g.strokeRect(L, top, w, l);
}

export function drawVehicle(g, size, stack, state, burning) {
  const st = stack.states[state], unit = size * .72 / stack.total_m * st.zoom, x = size / 2;
  const segs = st.segs.map((k) => stack.segments[k]), h = segs.reduce((a, q) => a + q.len, 0) * unit;
  let y = size / 2 + h / 2;                               // base of the stack; the sprite centre stays on the stack's centre
  const boost = stack.width_boost ?? WIDTH_BOOST, baseW = segs[0].w * unit * boost, own = OWN_BASE.has(segs[0].kind);
  g.clearRect(0, 0, size, size); g.lineJoin = 'round';
  if (st.boosters) {
    const b = st.boosters, bw = b.w * unit * boost, bl = b.len * unit, lw = Math.max(1.4, size * .007);
    for (const side of [-1, 1]) {
      const bx = x + side * (baseW / 2 + bw * .42);
      if (burning && b.burn) flame(g, bx, y + bw * .2, bw * .62, Math.min(size * .22, size - y - bw * .3), b.burn);
      g.fillStyle = shade(g, bx, bw * .8, COL.dark); g.strokeStyle = COL.outline; g.lineWidth = lw;
      g.beginPath(); g.moveTo(bx - bw * .2, y); g.lineTo(bx - bw * .34, y + bw * .22); g.lineTo(bx + bw * .34, y + bw * .22); g.lineTo(bx + bw * .2, y); g.closePath(); g.fill(); g.stroke();
      g.fillStyle = shade(g, bx, bw, COL.body); g.beginPath(); g.rect(bx - bw / 2, y - bl, bw, bl); g.fill(); g.stroke();
      g.beginPath(); g.moveTo(bx - bw / 2, y - bl); g.quadraticCurveTo(bx - bw * .42, y - bl - bw * .7, bx, y - bl - bw * .95); g.quadraticCurveTo(bx + bw * .42, y - bl - bw * .7, bx + bw / 2, y - bl); g.closePath(); g.fill(); g.stroke();
      g.fillStyle = COL.band; g.fillRect(bx - bw / 2 + 1, y - bl * .55, bw - 2, Math.max(2, bl * .05));
    }
  }
  g.strokeStyle = COL.outline; g.lineWidth = Math.max(1.4, size * .007);
  if (own) { if (burning && st.burn) flame(g, x, y, baseW * .8, size * .15, st.burn); }
  else {
    if (burning && st.burn) flame(g, x, y + baseW * .24, baseW * (st.burn === 'solid' ? .62 : .42), st.burn === 'solid' ? Math.min(size * .34, size - y - baseW * .3) : size * .13, st.burn);
    g.fillStyle = shade(g, x, baseW * .8, COL.dark);   // nozzle
    g.beginPath(); g.moveTo(x - baseW * .2, y); g.lineTo(x - baseW * .36, y + baseW * .28); g.lineTo(x + baseW * .36, y + baseW * .28); g.lineTo(x + baseW * .2, y); g.closePath(); g.fill(); g.stroke();
  }
  let prevW = baseW;
  for (const seg of segs) {
    const w = seg.w * unit * boost, l = seg.len * unit, top = y - l;
    if (w < prevW) {                                       // conical step down to a narrower stage
      const step = (prevW - w) * .55;
      g.fillStyle = shade(g, x, prevW, COL.body); g.beginPath(); g.moveTo(x - prevW / 2, y); g.lineTo(x - w / 2, y - step); g.lineTo(x + w / 2, y - step); g.lineTo(x + prevW / 2, y); g.closePath(); g.fill(); g.stroke();
      y -= step;
    }
    if (seg.kind === 'starship') { starship(g, x, y, w, l); break; }   // the ship is the top of its stack
    if (seg.kind === 'superheavy') { superHeavy(g, x, y, w, l); y = top; prevW = w; continue; }
    if (seg.kind === 'fairing') {
      const boat = (w - prevW) * 1.1, cyl = l * .46, nose = y - boat - cyl;
      g.fillStyle = shade(g, x, w, COL.body); g.beginPath(); g.moveTo(x - prevW / 2, y); g.lineTo(x - w / 2, y - boat); g.lineTo(x - w / 2, nose);
      g.quadraticCurveTo(x - w * .46, nose - l * .34, x, nose - l * .5 + boat * .2); g.quadraticCurveTo(x + w * .46, nose - l * .34, x + w / 2, nose); g.lineTo(x + w / 2, y - boat); g.lineTo(x + prevW / 2, y); g.closePath(); g.fill(); g.stroke();
      g.strokeStyle = 'rgba(40,50,60,.55)'; g.lineWidth = Math.max(1, size * .004); g.beginPath(); g.moveTo(x, y - boat); g.lineTo(x, nose - l * .22); g.stroke();
      break;
    }
    if (seg.kind === 'payload') {
      g.fillStyle = shade(g, x, w, COL.gold); g.strokeStyle = COL.outline; g.lineWidth = Math.max(1.4, size * .007);
      g.beginPath(); g.rect(x - w / 2, y - l, w, l); g.fill(); g.stroke();
      g.fillStyle = 'rgba(30,50,90,.85)'; g.fillRect(x - w * .9, y - l * .68, w * .34, l * .36); g.fillRect(x + w * .56, y - l * .68, w * .34, l * .36);   // stowed panels
      break;
    }
    const liquid = seg.kind === 'liquid';
    g.fillStyle = shade(g, x, w, liquid ? COL.dark : COL.body); g.strokeStyle = COL.outline; g.lineWidth = Math.max(1.4, size * .007);
    g.beginPath(); g.rect(x - w / 2, top, w, l); g.fill(); g.stroke();
    if (seg.band) { g.fillStyle = COL.band; g.fillRect(x - w / 2 + 1, y - l * .52, w - 2, Math.max(2, l * .055)); g.fillStyle = 'rgba(25,32,40,.92)'; g.fillRect(x - w / 2 + 1, y - l * .13, w - 2, Math.max(2, l * .05)); }
    if (!liquid) { g.fillStyle = 'rgba(25,32,40,.85)'; g.fillRect(x - w / 2 + 1, top, w - 2, Math.max(2, unit * .55)); }   // interstage ring
    y = top; prevW = w;
  }
}
