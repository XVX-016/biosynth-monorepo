export function encodeSharePayload(json: string): string {
  if (typeof window === 'undefined') return ''
  const encoded = window.btoa(unescape(encodeURIComponent(json)))
  return encoded.replace(/=+$/, '')
}

export function decodeSharePayload(encoded: string): string {
  if (typeof window === 'undefined') return ''
  const padded = encoded.padEnd(encoded.length + ((4 - (encoded.length % 4)) % 4), '=')
  return decodeURIComponent(escape(window.atob(padded)))
}



