import { OrbitControls } from "three/examples/jsm/controls/OrbitControls";

export class CameraController {
    controls: any;

    constructor(camera: any, domElement: any) {
        this.controls = new OrbitControls(camera, domElement);

        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.15;
        this.controls.zoomSpeed = 1.0;
    }

    update() {
        this.controls.update();
    }
}
