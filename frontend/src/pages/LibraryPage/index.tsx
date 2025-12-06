import LibrarySidebar from "./LibrarySidebar";
import LibraryContent from "./LibraryContent";
import { useLibrary } from "../../hooks/useLibrary";

export default function LibraryPage() {
    const library = useLibrary();

    return (
        <div className="flex h-full w-full">
            <LibrarySidebar upload={library.upload} />
            <LibraryContent
                molecules={library.molecules}
                loading={library.loading}
                remove={library.remove}
            />
        </div>
    );
}
