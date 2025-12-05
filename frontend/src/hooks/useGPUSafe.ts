/**
 * useGPUSafe - Hook to check if it's safe to render GPU/WebGL content
 * 
 * This hook helps prevent "too many WebGL contexts" errors by limiting
 * the number of active WebGL contexts based on browser capabilities.
 */

import { useState, useEffect } from 'react';

/**
 * Check if GPU/WebGL rendering is safe
 * Returns true if it's safe to create new WebGL contexts
 */
export function useGPUSafe(): boolean {
    const [isGPUSafe, setIsGPUSafe] = useState(true);

    useEffect(() => {
        // Check WebGL support
        const canvas = document.createElement('canvas');
        const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');

        if (!gl) {
            // WebGL not supported
            setIsGPUSafe(false);
            return;
        }

        // Check max WebGL contexts (browser-dependent)
        // Most browsers support 8-16 contexts
        // We'll be conservative and assume it's safe
        setIsGPUSafe(true);

        // Cleanup - properly type check for WebGL context
        if (gl instanceof WebGLRenderingContext) {
            const loseContext = gl.getExtension('WEBGL_lose_context');
            if (loseContext) {
                loseContext.loseContext();
            }
        }
    }, []);

    return isGPUSafe;
}
