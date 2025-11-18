import React from 'react';
import clsx from 'clsx';

export type ButtonProps = React.ButtonHTMLAttributes<HTMLButtonElement> & {
	variant?: 'primary' | 'secondary' | 'danger';
	size?: 'sm' | 'md' | 'lg';
};

export default function Button({
	variant = 'primary',
	size = 'md',
	className,
	disabled,
	children,
	...rest
}: ButtonProps) {
	const variantClasses = {
		primary: 'bg-black text-white border border-black hover:bg-darkGrey hover:shadow-neon disabled:opacity-50 font-semibold interactive-glow',
		secondary: 'bg-white text-black border border-lightGrey hover:bg-offwhite hover:shadow-neon disabled:opacity-50 interactive-glow',
		danger: 'bg-white text-black border border-lightGrey hover:bg-offwhite hover:shadow-neon disabled:opacity-50 interactive-glow',
	};

	const sizeClasses = {
		sm: 'px-2 py-1 text-sm',
		md: 'px-4 py-2',
		lg: 'px-6 py-3 text-lg',
	};

	return (
		<button
			{...rest}
			disabled={disabled}
			className={clsx(
				'inline-flex items-center justify-center rounded-lg font-medium transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-darkGrey/20',
				variantClasses[variant],
				sizeClasses[size],
				className
			)}
		>
			{children}
		</button>
	);
}


