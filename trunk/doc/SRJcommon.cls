%  Common classes and definitions for SRJhomework etc.
%  by Seth R. Johnson
% licensed under Creative Commons (attribution): http://creativecommons.org/licenses/by/3.0/us/ 

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{SRJcommon}
[2008/09/02 v1.00 Common classes and definitions, Seth R. Johnson]

%%% OPTIONS %%%

\newif\ifusesubfig
\usesubfigfalse

\newif\ifusetinyunits
\usetinyunitsfalse

\newif\ifusefancyvrb
\usefancyvrbfalse

\newif\ifusecolor
\usecolorfalse

\DeclareOption{tinyunits}{%
	\usetinyunitstrue
	\usecolortrue}
\DeclareOption{fancyvrb}{%
	\usefancyvrbtrue}
\DeclareOption{subfig}{%
	\usesubfigtrue}
\DeclareOption{color}{%
	\usecolortrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

%%% REQUIRED CLASSES %%%
\LoadClass{article}
\RequirePackage{geometry}
\RequirePackage{fancyhdr}

\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{amsthm}
\RequirePackage{bm}      %bold math symbols for vectors
\RequirePackage[pdftex]{graphicx}
\RequirePackage{caption}

\captionsetup{width=.75\textwidth,font=small,labelfont=bf}
%\captionsetup{width=.75\textwidth,font=small,format=hang}

%%% CLASS NAME AND TITLES %%%
\def\@classname{\PackageWarning{No `classname' given.}}
\newcommand{\classname}[1]{\def\@classname{#1}}

%%% MISC PACKAGE DEFS
% colors
\ifusecolor
  \RequirePackage{color}
  \definecolor{lightgray}{gray}{0.5}
  \definecolor{lightblue}{rgb}{0.4,0.4,0.80} 
  \definecolor{red}{rgb}{1,0,0} 
\fi

% set up subfigures
\ifusesubfig
	%\typeout{Using subfig}
	\RequirePackage{subfigure}%
%
	\renewcommand{\thesubfigure}{\alph{subfigure}}%
	\renewcommand{\@thesubfigure}{{\subcaplabelfont(\thesubfigure)}\space}%
\fi

%set up fancyvrb
\ifusefancyvrb
	\usepackage{fancyvrb}%
%
	\newcommand{\VerbCode}[2]{\fvset{frame=single,%
		fontsize=\small,%
		xleftmargin=0.5in,%
		xrightmargin=0.5in,%
		label={#1},%
		labelposition=topline%
	}%
	\VerbatimInput{#2}%
	}%
\fi

% REQUIRED CLASSES (load hyperref last because of dependencies)
\RequirePackage[unicode,pdftex]{hyperref}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        COMMAND SHORTCUTS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcuts in regular text
\newcommand{\degs}{\ensuremath{^\circ}}
\newcommand{\EE}[1]{\ensuremath{\times 10^{#1}}}
\newcommand{\super}[1]{\ensuremath{{}^{\text{#1}}}}
\newcommand{\ttimes}{\ensuremath{{}\times{}}}
\newcommand{\SN}{S$_\text{N}$}

% Program names
\newcommand{\cpp}{\textsf{C\raisebox{0.2ex}{++}}}

% General math shortcuts
\newcommand{\ud}{\, \mathrm{d}}
\newcommand{\Ud}{\mathrm{d}}
\newcommand{\pder}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\oder}[2]{\frac{\mathrm{d} #1}{\mathrm{d} #2}}
\newcommand{\lra}{ \quad \Longrightarrow \quad }
\newcommand{\expp}[1]{\exp\!\left( {#1} \right)}

% Numerical Linear Algebra
\newcommand{\conj}{^{\ast}}
\newcommand{\norm}[1]{\| #1 \|}
\newcommand{\abs}[1]{| #1 |}
%\newcommand{\norm}[1]{\left\| #1 \right\|}
%\newcommand{\abs}[1]{\left| #1 \right|}
\newcommand{\eps}{\varepsilon}
\DeclareMathOperator{\fl}{fl}

\DeclareMathOperator{\acosh}{arccosh} 

% Define a command to write a nice-looking element, e.g. 4,2 He
\newcommand{\elem}[3]{\ensuremath{{}^{{#1}}_{{#2}}\mathrm{{#3}}}}

%Units: first is the tiny blue offset label, second is the more formal way
\ifusetinyunits%
	\newcommand{\units}[1] {\raisebox{0.3ex}{\color{lightblue}\tiny{[{#1}]}}}%
	\newcommand{\ohms}{\units{ohm}}%
\else%
	\newcommand{\units}[1] {\:\text{#1}}%
%	\newcommand{\ohms}{\ensuremath{\Omega}}%
	\newcommand{\ohms}{\units{$\Omega$}}%
\fi

% Vector definitions
\renewcommand{\vec}[1]{\ensuremath{\bm{#1}}}
%\renewcommand{\vec}[1]{\ensuremath{\boldsymbol{#1}}} %use italic bold font
%\renewcommand{\vec}[1]{\ensuremath{\mathbf{#1}}} %use bold math
%\newcommand{\vd}{\bullet}
%\newcommand{\vd}{\cdot}
\newcommand{\vd}{\vec{\cdot}}
\newcommand{\del}{\vec{\nabla}}
\newcommand{\op}[1]{\ensuremath{\mathsf{#1}}}

%\newcommand{\abr}[1]{\langle {#1} \rangle}
\newcommand{\abr}[1]{\left\langle {#1} \right\rangle}

% Shielding
\newcommand{\resp}{\mathcal{R}} %dose response

%alternate definition adds a bit more room between cases
\newenvironment{Cases}{%
  \left\{\def\arraystretch{1.5}%
  \array{@{}l@{\quad}l@{}}%
}{%
  \endarray\right.%
} 

%% topbox is useful in extended definitions of math terms inside an align
\newcommand{\topbox}[1]{\parbox[t]{0.6\columnwidth}{\raggedright{}#1}}

% commands to make text in math mode appear as zero-width (better-looking integrals/sums, e.g.)
% from mathmode.pdf page 74
\def\mathllap{\mathpalette\mathllapinternal}
	\def\mathllapinternal#1#2{%
	\llap{$\mathsurround=0pt#1{#2}$}%
}
\def\clap#1{\hbox to 0pt{\hss#1\hss}}%
\def\mathclap{\mathpalette\mathclapinternal}%
\def\mathclapinternal#1#2{%
	\clap{$\mathsurround=0pt#1{#2}$}%
}
\def\mathrlap{\mathpalette\mathrlapinternal}%
\def\mathrlapinternal#1#2{%
	\rlap{$\mathsurround=0pt#1{#2}$}%
}

