import React from 'react';
import clsx from 'clsx';

export type CardProps = {
	header?: React.ReactNode;
	footer?: React.ReactNode;
	children: React.ReactNode;
	className?: string;
};

export default function Card({ header, footer, children, className }: CardProps) {
	return (
		<section
			className={clsx(
				'frosted-glass rounded-xl shadow-glass border border-chrome/20',
				className
			)}
		>
			{header && <div className="px-5 py-4 border-b border-chrome/20">{header}</div>}
			<div className="px-5 py-4 text-ivory">{children}</div>
			{footer && <div className="px-5 py-4 border-t border-chrome/20">{footer}</div>}
		</section>
	);
}


