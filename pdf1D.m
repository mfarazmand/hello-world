function [xx,pdf_f]=pdf1D(f,N)

[Ni,xx]=hist(f,N);

dx=xx(2)-xx(1);
pdf_f = Ni/(sum(Ni)*dx);