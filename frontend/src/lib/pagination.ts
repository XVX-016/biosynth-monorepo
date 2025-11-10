export function paginate<T>(items: T[], page: number, pageSize: number): { data: T[]; totalPages: number } {
	const totalPages = Math.max(1, Math.ceil((items?.length ?? 0) / pageSize));
	const current = Math.min(Math.max(1, page), totalPages);
	const start = (current - 1) * pageSize;
	const end = start + pageSize;
	return { data: (items ?? []).slice(start, end), totalPages };
}


