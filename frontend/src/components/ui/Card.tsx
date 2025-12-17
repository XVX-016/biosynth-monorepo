import React from 'react';
import clsx from 'clsx';

export type CardProps = {
	header?: React.ReactNode;
	footer?: React.ReactNode;
	children: React.ReactNode;
	className?: string;
	style?: React.CSSProperties;
};

export default function Card({ header, footer, children, className, style }: CardProps) {
	return (
		<section
			className={clsx(
				'bg-white rounded-xl shadow-neon border border-lightGrey',
				className
			)}
			style={style}
		>
			{header && <div className="px-5 py-4 border-b border-lightGrey">{header}</div>}
			<div className="px-5 py-4 text-black">{children}</div>
			{footer && <div className="px-5 py-4 border-t border-lightGrey">{footer}</div>}
		</section>
	);
}


