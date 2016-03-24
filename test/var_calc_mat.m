f=rand(10,1);
var(f)
sum(f.^2)/length(f)-(sum(f)/length(f))^2
%ones(1,10)*(diag(f.^2)/9-f*f'/10)*ones(10,1)
(f'*f-ones(1,10)*(f*f')*ones(10,1)/length(f))/length(f)
ones(1,10)*(diag(f.^2)-(f*f')/length(f))/length(f)*ones(10,1)
%F=(diag(f.^2)-(f*f')/length(f))
[B,C]=eigs(f*f'/length(f))

