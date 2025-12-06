export const CH4 = {
    name: "Methane",
    atoms: [
        { id: "C1", element: "C", x: 0, y: 0, z: 0 },
        { id: "H1", element: "H", x: 1, y: 1, z: 1 },
        { id: "H2", element: "H", x: -1, y: -1, z: 1 },
        { id: "H3", element: "H", x: 1, y: -1, z: -1 },
        { id: "H4", element: "H", x: -1, y: 1, z: -1 }
    ],
    bonds: [
        { a: "C1", b: "H1", order: 1 },
        { a: "C1", b: "H2", order: 1 },
        { a: "C1", b: "H3", order: 1 },
        { a: "C1", b: "H4", order: 1 }
    ]
};

export const BENZENE = {
    name: "Benzene",
    atoms: [
        { id: "C1", element: "C", x: 1.396, y: 0, z: 0 },
        { id: "C2", element: "C", x: 0.698, y: 1.209, z: 0 },
        { id: "C3", element: "C", x: -0.698, y: 1.209, z: 0 },
        { id: "C4", element: "C", x: -1.396, y: 0, z: 0 },
        { id: "C5", element: "C", x: -0.698, y: -1.209, z: 0 },
        { id: "C6", element: "C", x: 0.698, y: -1.209, z: 0 }
    ],
    bonds: [
        { a: "C1", b: "C2", order: 1 },
        { a: "C2", b: "C3", order: 2 },
        { a: "C3", b: "C4", order: 1 },
        { a: "C4", b: "C5", order: 2 },
        { a: "C5", b: "C6", order: 1 },
        { a: "C6", b: "C1", order: 2 }
    ]
};
