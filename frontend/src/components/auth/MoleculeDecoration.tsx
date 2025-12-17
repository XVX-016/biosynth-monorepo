/**
 * MoleculeDecoration - Subtle animated SVG molecule
 * 
 * Lightweight, GPU-accelerated SVG animation
 * No React re-renders, no JS timers
 * Hidden on mobile for performance
 */
export default function MoleculeDecoration() {
  return (
    <svg
      viewBox="0 0 300 300"
      className="w-full h-full opacity-[0.08]"
      aria-hidden="true"
      xmlns="http://www.w3.org/2000/svg"
    >
      {/* Central atom */}
      <circle cx="150" cy="150" r="12" fill="currentColor" className="text-black" />
      
      {/* Outer atoms */}
      <circle cx="150" cy="80" r="8" fill="currentColor" className="text-black" />
      <circle cx="80" cy="150" r="8" fill="currentColor" className="text-black" />
      <circle cx="220" cy="150" r="8" fill="currentColor" className="text-black" />
      <circle cx="150" cy="220" r="8" fill="currentColor" className="text-black" />
      <circle cx="105" cy="105" r="7" fill="currentColor" className="text-black" />
      <circle cx="195" cy="105" r="7" fill="currentColor" className="text-black" />
      <circle cx="195" cy="195" r="7" fill="currentColor" className="text-black" />
      <circle cx="105" cy="195" r="7" fill="currentColor" className="text-black" />

      {/* Bonds */}
      <line x1="150" y1="158" x2="150" y2="88" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="158" y1="150" x2="228" y2="150" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="142" y1="150" x2="72" y2="150" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="150" y1="142" x2="150" y2="212" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="157" y1="157" x2="112" y2="112" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="143" y1="157" x2="188" y2="112" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="157" y1="143" x2="112" y2="188" stroke="currentColor" strokeWidth="2" className="text-black" />
      <line x1="143" y1="143" x2="188" y2="188" stroke="currentColor" strokeWidth="2" className="text-black" />

      {/* Gentle bounce/float animation */}
      <animateTransform
        attributeName="transform"
        type="translate"
        values="0,0; 0,-10; 0,0; 0,5; 0,0"
        dur="8s"
        repeatCount="indefinite"
      />
      {/* Slow rotation */}
      <animateTransform
        attributeName="transform"
        type="rotate"
        from="0 150 150"
        to="360 150 150"
        dur="25s"
        repeatCount="indefinite"
        additive="sum"
      />
    </svg>
  );
}

