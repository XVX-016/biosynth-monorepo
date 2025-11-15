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
		primary: 'bg-plasma-neon text-ionBlack hover:shadow-neon-sm disabled:opacity-50 font-semibold',
		secondary: 'bg-frostedGlass text-ivory border border-chrome/30 hover:border-neonCyan/50 hover:shadow-neon-sm disabled:opacity-50 backdrop-blur-sm',
		danger: 'bg-violetEdge/20 text-violetEdge border border-violetEdge/30 hover:bg-violetEdge/30 disabled:opacity-50',
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
				'inline-flex items-center justify-center rounded-lg font-medium transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-neonCyan/50',
				variantClasses[variant],
				sizeClasses[size],
				className
			)}
		>
			{children}
		</button>
	);
}


