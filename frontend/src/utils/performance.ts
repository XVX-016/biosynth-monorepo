/**
 * Performance monitoring utilities
 * 
 * Phase 16: Full Regression Test & Final Cleanup
 * 
 * Provides performance profiling and monitoring.
 */

/**
 * Measure execution time of a function
 */
export function measurePerformance<T>(
  name: string,
  fn: () => T,
  log: boolean = true
): T {
  const start = performance.now()
  const result = fn()
  const end = performance.now()
  const duration = end - start

  if (log) {
    console.log(`[Performance] ${name}: ${duration.toFixed(2)}ms`)
  }

  return result
}

/**
 * Measure async execution time
 */
export async function measureAsyncPerformance<T>(
  name: string,
  fn: () => Promise<T>,
  log: boolean = true
): Promise<T> {
  const start = performance.now()
  const result = await fn()
  const end = performance.now()
  const duration = end - start

  if (log) {
    console.log(`[Performance] ${name}: ${duration.toFixed(2)}ms`)
  }

  return result
}

/**
 * Monitor frame rate
 */
export class FrameRateMonitor {
  private frames: number[] = []
  private startTime: number = 0
  private frameCount: number = 0
  private rafId: number | null = null
  private callback: ((fps: number) => void) | null = null

  start(callback?: (fps: number) => void) {
    this.callback = callback || null
    this.frames = []
    this.frameCount = 0
    this.startTime = performance.now()
    this.measure()
  }

  private measure = () => {
    const now = performance.now()
    this.frameCount++

    if (now >= this.startTime + 1000) {
      const fps = this.frameCount
      this.frames.push(fps)
      if (this.frames.length > 60) {
        this.frames.shift()
      }

      if (this.callback) {
        this.callback(fps)
      }

      this.frameCount = 0
      this.startTime = now
    }

    this.rafId = requestAnimationFrame(this.measure)
  }

  stop() {
    if (this.rafId !== null) {
      cancelAnimationFrame(this.rafId)
      this.rafId = null
    }
  }

  getAverageFPS(): number {
    if (this.frames.length === 0) return 0
    return this.frames.reduce((a, b) => a + b, 0) / this.frames.length
  }

  getMinFPS(): number {
    if (this.frames.length === 0) return 0
    return Math.min(...this.frames)
  }

  getMaxFPS(): number {
    if (this.frames.length === 0) return 0
    return Math.max(...this.frames)
  }
}

/**
 * Monitor memory usage (if available)
 */
export function getMemoryUsage(): {
  used: number
  total: number
  percentage: number
} | null {
  if ('memory' in performance) {
    const memory = (performance as any).memory
    return {
      used: memory.usedJSHeapSize,
      total: memory.totalJSHeapSize,
      percentage: (memory.usedJSHeapSize / memory.totalJSHeapSize) * 100,
    }
  }
  return null
}

/**
 * Check WebGL context capabilities
 */
export function checkWebGLCapabilities(): {
  supported: boolean
  maxTextureSize: number
  maxVertexAttribs: number
  maxVaryingVectors: number
  maxFragmentUniformVectors: number
  maxVertexUniformVectors: number
} {
  const canvas = document.createElement('canvas')
  const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl')

  if (!gl) {
    return {
      supported: false,
      maxTextureSize: 0,
      maxVertexAttribs: 0,
      maxVaryingVectors: 0,
      maxFragmentUniformVectors: 0,
      maxVertexUniformVectors: 0,
    }
  }

  return {
    supported: true,
    maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE),
    maxVertexAttribs: gl.getParameter(gl.MAX_VERTEX_ATTRIBS),
    maxVaryingVectors: gl.getParameter(gl.MAX_VARYING_VECTORS),
    maxFragmentUniformVectors: gl.getParameter(gl.MAX_FRAGMENT_UNIFORM_VECTORS),
    maxVertexUniformVectors: gl.getParameter(gl.MAX_VERTEX_UNIFORM_VECTORS),
  }
}

