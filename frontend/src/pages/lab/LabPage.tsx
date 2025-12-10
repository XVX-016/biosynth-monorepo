import React, { useEffect } from 'react';
import LabLayout from '../../components/lab/LabLayout';
import { EditorProvider, useEditor } from '../../context/EditorContext';

function LabLoader() {
    const { dispatch } = useEditor();
    useEffect(() => {
        const params = new URLSearchParams(window.location.search);
        const molId = params.get('molId');
        if (!molId) return;
        (async () => {
            try {
                const res = await fetch(`/api/library/molecule/${molId}`);
                if (!res.ok) throw new Error('Not found');
                const json = await res.json();
                // expected: { atoms: [...], bonds: [...] }
                dispatch({ type: 'LOAD_MOLECULE', payload: { atoms: json.atoms, bonds: json.bonds } });
            } catch (e) {
                console.warn('Failed to load molecule', e);
            }
        })();
    }, [dispatch]);
    return null;
}

export default function LabPage() {
    return (
        <EditorProvider>
            <LabLoader />
            <LabLayout />
        </EditorProvider>
    );
}
