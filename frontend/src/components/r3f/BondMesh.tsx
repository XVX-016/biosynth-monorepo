import React from 'react'
import * as THREE from 'three'

interface BondMeshProps {
  from: [number, number, number]
  to: [number, number, number]
}

export default function BondMesh({ from, to }: BondMeshProps) {
  const vFrom = new THREE.Vector3(...from)
  const vTo = new THREE.Vector3(...to)
  const diff = new THREE.Vector3().subVectors(vTo, vFrom)
  const length = diff.length()
  const mid = new THREE.Vector3().addVectors(vFrom, vTo).multiplyScalar(0.5)
  const q = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0), 
    diff.clone().normalize()
  )
  
  return (
    <mesh 
      position={[mid.x, mid.y, mid.z]} 
      quaternion={[q.x, q.y, q.z, q.w]}
    >
      <cylinderGeometry args={[0.14, 0.14, length, 16]} />
      <meshStandardMaterial 
        color={0x8f9aa3} 
        metalness={0.8} 
        roughness={0.4}
      />
    </mesh>
  )
}

