import { create } from 'zustand';

type ProfileState = {
	name: string;
	avatarUrl: string | null;
	setName: (name: string) => void;
	setAvatarUrl: (url: string | null) => void;
	loadFromStorage: () => void;
};

const STORAGE_KEY = 'biosynth.profile';

export const useProfileStore = create<ProfileState>((set, get) => ({
	name: 'Researcher',
	avatarUrl: null,
	setName: (name: string) => {
		set({ name });
		try {
			const current = get();
			localStorage.setItem(
				STORAGE_KEY,
				JSON.stringify({ name, avatarUrl: current.avatarUrl })
			);
		} catch {
			// ignore
		}
	},
	setAvatarUrl: (url: string | null) => {
		set({ avatarUrl: url });
		try {
			const current = get();
			localStorage.setItem(
				STORAGE_KEY,
				JSON.stringify({ name: current.name, avatarUrl: url })
			);
		} catch {
			// ignore
		}
	},
	loadFromStorage: () => {
		try {
			const raw = localStorage.getItem(STORAGE_KEY);
			if (!raw) return;
			const parsed = JSON.parse(raw) as { name?: string; avatarUrl?: string | null };
			set({
				name: parsed.name ?? 'Researcher',
				avatarUrl: parsed.avatarUrl ?? null,
			});
		} catch {
			// ignore
		}
	},
}));


