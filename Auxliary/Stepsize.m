function lamb=Stepsize(p1,p2,p3,p4)
% if p1<=1e-40
%     lamb=1;
% else

if p2<0 && p3>0 &&(p2*p3<=p1*p4)
    lamb=min([1,-p4/p3]);
elseif p1+p2+p3+p4<=0
    lamb=1;
else
    lamb=-p4/(p1+p2+p3);
end
%     elseif p1+p2+p3+p4>0
%         lamb=-p4/(p1+p2+p3);
%     else
%         lamb=1;
%         % elseif p1+p2+p3+p4<=0
%         %     lamb=1;
%         % else
%     end
%
% lamb=min([lamb,1]);

% if p2>=0 || p2*p3>=0 ||(p2*p3<0&&3*p1*p4<=p2*p3)
%     if p1+p2+p3+p4<=0
%         lamb=1;
%     else
%         lamb=-p4/(p1+p2+p3);
%     end
% else
%     lamb=-p4/p3;
% end