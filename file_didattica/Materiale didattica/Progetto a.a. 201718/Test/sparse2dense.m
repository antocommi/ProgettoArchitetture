function P = sparse2dense(G)
    c = 0.85;
    n = max(max(G));
    deg = hist(G(:,1),1:n)';
    P = spconvert([G 1./deg(G(:,1))]);
    if (diff(size(P))~=0)
        P(n,n) = 0;
    end
    P(deg==0,:) = 1/n;
    P = full(P);
    P = c*P+(1-c)/n;
end
