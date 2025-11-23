/**
 * GPU Safety Hook
 * Prevents heavy 3D rendering on mobile devices and low-end hardware
 */
import { useState, useEffect } from 'react';

export function useGPUSafe(): boolean {
  const [isGPUSafe, setIsGPUSafe] = useState(false);

  useEffect(() => {
    // Check if device is desktop-sized (min-width: 900px)
    const mediaQuery = window.matchMedia('(min-width: 900px)');
    
    const updateGPUSafe = () => {
      setIsGPUSafe(mediaQuery.matches);
    };

    // Set initial value
    updateGPUSafe();

    // Listen for changes
    mediaQuery.addEventListener('change', updateGPUSafe);

    return () => {
      mediaQuery.removeEventListener('change', updateGPUSafe);
    };
  }, []);

  return isGPUSafe;
}

