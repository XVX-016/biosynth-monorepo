import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';

export default class Renderer3D {
    container: HTMLElement;
    renderer: THREE.WebGLRenderer;
    scene: THREE.Scene;
    camera: THREE.PerspectiveCamera;
    controls: any;
    raf = 0;

    // track nodes by id
    atomsMap = new Map<string, THREE.Object3D>();
    bondsMap = new Map<string, THREE.Mesh>();

    constructor(container: HTMLElement) {
        this.container = container;
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        container.appendChild(this.renderer.domElement);

        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0b0b0d);

        this.camera = new THREE.PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 1000);
        this.camera.position.set(0, 0, 50);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;

        const am = new THREE.AmbientLight(0xffffff, 0.6);
        const dir = new THREE.DirectionalLight(0xffffff, 0.8);
        dir.position.set(10, 10, 10);
        this.scene.add(am, dir);

        const g = new THREE.GridHelper(200, 50, 0x222222, 0x111111);
        g.rotation.x = Math.PI / 2;
        g.position.y = -20;
        this.scene.add(g);

        window.addEventListener('resize', this.onWindowResize);
    }

    start = () => {
        const loop = (t = 0) => {
            this.controls.update();
            this.renderer.render(this.scene, this.camera);
            this.raf = requestAnimationFrame(loop);
        };
        loop();
    }

    onWindowResize = () => {
        const w = this.container.clientWidth;
        const h = this.container.clientHeight;
        this.camera.aspect = w / h;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(w, h);
    }

    addAtomSphere(id: string, element: string, pos: [number, number, number]) {
        const geo = new THREE.SphereGeometry(0.8, 16, 16);
        const mat = new THREE.MeshStandardMaterial({ metalness: 0.1, roughness: 0.6 });
        const m = new THREE.Mesh(geo, mat);
        m.position.set(...pos);
        (m as any).userData = { id, type: 'atom', element };
        this.scene.add(m);
        this.atomsMap.set(id, m);
        return m;
    }

    // create bond as a cylinder between two atom positions
    addBondCylinder(id: string, aPos: [number, number, number], bPos: [number, number, number], order = 1) {
        const a = new THREE.Vector3(...aPos);
        const b = new THREE.Vector3(...bPos);
        const mid = new THREE.Vector3().addVectors(a, b).multiplyScalar(0.5);
        const dir = new THREE.Vector3().subVectors(b, a);
        const len = dir.length();

        const radius = 0.18 * order; // visual thickness
        const geo = new THREE.CylinderGeometry(radius, radius, len, 12, 1);
        const mat = new THREE.MeshStandardMaterial({ metalness: 0.2, roughness: 0.4 });
        const mesh = new THREE.Mesh(geo, mat);

        // orient the cylinder
        mesh.position.copy(mid);
        mesh.lookAt(b);
        mesh.rotateX(Math.PI / 2);

        (mesh as any).userData = { id, type: 'bond', order };
        this.scene.add(mesh);
        this.bondsMap.set(id, mesh);
        return mesh;
    }

    // animated creation: scale from 0 -> 1 on Y
    async addBondAnimated(id: string, aPos: [number, number, number], bPos: [number, number, number], order = 1, duration = 220) {
        const mesh = this.addBondCylinder(id, aPos, bPos, order);
        mesh.scale.set(1, 0.0001, 1);
        const start = performance.now();
        return new Promise<void>((resolve) => {
            const step = (now: number) => {
                const p = Math.min(1, (now - start) / duration);
                mesh.scale.y = p;
                if (p < 1) requestAnimationFrame(step);
                else resolve();
            };
            requestAnimationFrame(step);
        });
    }

    removeBond(id: string) {
        const m = this.bondsMap.get(id);
        if (!m) return;
        this.scene.remove(m);
        this.bondsMap.delete(id);
    }

    onPointerDown(e: PointerEvent) {
        const rect = this.renderer.domElement.getBoundingClientRect();
        const x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
        const y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
        const ray = new THREE.Raycaster();
        ray.setFromCamera(new THREE.Vector2(x, y), this.camera);
        const hits = ray.intersectObjects(this.scene.children, true);
        if (hits.length) {
            const obj = hits[0].object;
            console.log('clicked', (obj as any).userData);
        }
    }

    dispose() {
        cancelAnimationFrame(this.raf);
        this.renderer.dispose();
        window.removeEventListener('resize', this.onWindowResize);
        this.container.removeChild(this.renderer.domElement);
    }
}
