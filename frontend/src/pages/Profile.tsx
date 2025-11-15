import React, { useEffect, useMemo, useState } from 'react';
import { listMolecules } from '../lib/api';
import { useProfileStore } from '../store/profileStore';

export default function Profile() {
	const name = useProfileStore((s) => s.name);
	const setName = useProfileStore((s) => s.setName);
	const load = useProfileStore((s) => s.loadFromStorage);
	const [count, setCount] = useState(0);
	useEffect(() => {
		load();
		let cancelled = false;
		(async () => {
			try {
				const res = await listMolecules(100);
				if (!cancelled) setCount(res?.length ?? 0);
			} catch {
				// ignore for now
			}
		})();
		return () => {
			cancelled = true;
		};
	}, []);
	const recentCount = useMemo(() => Math.min(count, 5), [count]);

	return (
		<div className="grid grid-cols-12 gap-4">
			<div className="col-span-12 lg:col-span-4">
				<div className="frosted-glass rounded-xl shadow-glass border border-chrome/20 p-6">
					<div className="flex items-center gap-4">
						<div className="w-16 h-16 rounded-full bg-spaceGrey border border-chrome/30" />
						<div>
							<h2 className="text-xl font-semibold text-ivory">Your Profile</h2>
							<p className="text-chrome">Local profile</p>
						</div>
					</div>
					<div className="mt-6">
						<label className="block text-sm text-chrome mb-1">Display Name</label>
						<input
							type="text"
							placeholder="Your name"
							value={name}
							onChange={(e) => setName(e.target.value)}
							className="w-full rounded-lg border border-chrome/20 bg-frostedGlass text-ivory px-3 py-2 outline-none focus:ring-2 focus:ring-neonCyan/50 focus:border-neonCyan/50"
						/>
					</div>
				</div>
			</div>
			<div className="col-span-12 lg:col-span-8">
				<div className="frosted-glass rounded-xl shadow-glass border border-chrome/20 p-6">
					<h3 className="text-lg font-semibold text-ivory">Summary</h3>
					<div className="mt-4 grid grid-cols-2 gap-4">
						<div className="rounded-xl border border-chrome/20 frosted-glass p-4">
							<div className="text-sm text-chrome">Saved molecules</div>
							<div className="text-2xl font-bold text-ivory">{count}</div>
						</div>
						<div className="rounded-xl border border-chrome/20 frosted-glass p-4">
							<div className="text-sm text-chrome">Recent activity</div>
							<div className="text-2xl font-bold text-ivory">{recentCount}</div>
						</div>
					</div>
				</div>
			</div>
		</div>
	);
}


