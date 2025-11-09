import React, { Suspense } from 'react'
import { Canvas } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import AtomMesh from './r3f/AtomMesh'
import BondMesh from './r3f/BondMesh'

export default function MoleculeViewer(){
  return (
    <div style={{height:420}} className="rounded-lg overflow-hidden bg-aluminum-light">
      <Canvas dpr={[1,2]} camera={{ position: [0,0,12], fov:45 }}>
        <ambientLight intensity={0.6}/>
        <directionalLight position={[5,10,7]} intensity={0.6} />
        <Suspense fallback={null}>
          {/* Simple methane demo - CH4 */}
          <AtomMesh position={[0,0,0]} element="C"/>
          <AtomMesh position={[2.4,1.2,0]} element="H"/>
          <AtomMesh position={[-2.4,1.2,0]} element="H"/>
          <AtomMesh position={[0,-1.4,2.0]} element="H"/>
          <AtomMesh position={[0,-1.4,-2.0]} element="H"/>
          <BondMesh from={[0,0,0]} to={[2.4,1.2,0]} />
          <BondMesh from={[0,0,0]} to={[-2.4,1.2,0]} />
          <BondMesh from={[0,0,0]} to={[0,-1.4,2.0]} />
          <BondMesh from={[0,0,0]} to={[0,-1.4,-2.0]} />
        </Suspense>
        <OrbitControls enablePan={true} enableZoom={true} enableRotate={true} />
      </Canvas>
    </div>
  )
}

