import React from 'react';
import clsx from 'clsx';

export type PanelProps = {
	children: React.ReactNode;
	className?: string;
};

export default function Panel({ children, className }: PanelProps) {
	return (
		<div
			className={clsx(
				'frosted-glass rounded-xl shadow-glass border border-chrome/20 p-4 text-ivory',
				className
			)}
		>
			{children}
		</div>
	);
}


