import React from "react";
import OptimizeButton from "./OptimizeButton";
import PredictButton from "./PredictButton";
import AutoBondButton from "./tools/AutoBondButton";
import AddAtomButton from "./tools/AddAtomButton";
import AddBondButton from "./tools/AddBondButton";
import ExportButton from "./ExportButton";

export default function Toolbar() {
    return (
        <div className="flex gap-2 p-2 bg-white border-b flex-wrap">
            <AddAtomButton />
            <AddBondButton />
            <AutoBondButton />
            <OptimizeButton />
            <PredictButton />
            <ExportButton />
        </div>
    );
}
