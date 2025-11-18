import React, { useEffect, useMemo, useState } from 'react';
import { listMolecules } from '../lib/api';
import { useProfileStore } from '../store/profileStore';
import Card from '../components/ui/Card';

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
		<div className="grid grid-cols-12 gap-6">
			<div className="col-span-12 lg:col-span-4">
				<Card className="p-6">
					<div className="flex items-center gap-4">
						<div className="w-16 h-16 rounded-full bg-lightGrey border border-lightGrey" />
						<div>
							<h2 className="text-xl font-semibold text-black">Your Profile</h2>
							<p className="text-midGrey">Local profile</p>
						</div>
					</div>
					<div className="mt-6">
						<label className="block text-sm text-darkGrey mb-1">Display Name</label>
						<input
							type="text"
							placeholder="Your name"
							value={name}
							onChange={(e) => setName(e.target.value)}
							className="w-full rounded-lg border border-lightGrey bg-white text-black px-3 py-2 outline-none focus:ring-2 focus:ring-darkGrey/20 focus:border-darkGrey"
						/>
					</div>
				</Card>
			</div>
			<div className="col-span-12 lg:col-span-8">
				<Card className="p-6">
					<h3 className="text-lg font-semibold text-black">Summary</h3>
					<div className="mt-4 grid grid-cols-2 gap-4">
						<Card className="p-4">
							<div className="text-sm text-midGrey">Saved molecules</div>
							<div className="text-2xl font-bold text-black">{count}</div>
						</Card>
						<Card className="p-4">
							<div className="text-sm text-midGrey">Recent activity</div>
							<div className="text-2xl font-bold text-black">{recentCount}</div>
						</Card>
					</div>
				</Card>
			</div>
		</div>
	);
}


