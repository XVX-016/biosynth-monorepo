import MoleculeSearch from "./MoleculeSearch";
import MoleculeFilters from "./MoleculeFilters";
import MoleculeUploader from "./MoleculeUploader";

interface LibrarySidebarProps {
    upload: (mol: any) => Promise<void>;
}

export default function LibrarySidebar({ upload }: LibrarySidebarProps) {
    return (
        <div className="w-72 bg-white border-r p-4 flex flex-col gap-4">
            <MoleculeSearch />
            <MoleculeFilters />
            <MoleculeUploader onUpload={upload} />
        </div>
    );
}
