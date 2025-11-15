import React from 'react';
import clsx from 'clsx';

export type IconButtonProps = React.ButtonHTMLAttributes<HTMLButtonElement> & {
	'aria-label': string;
};

export default function IconButton({ className, children, ...rest }: IconButtonProps) {
	return (
		<button
			type="button"
			{...rest}
			className={clsx(
				'inline-flex h-9 w-9 items-center justify-center rounded-lg border border-chrome/20 bg-frostedGlass text-chrome hover:text-ivory hover:border-neonCyan/30 focus:outline-none focus:ring-2 focus:ring-neonCyan/50 transition-all',
				className
			)}
		>
			{children}
		</button>
	);
}


