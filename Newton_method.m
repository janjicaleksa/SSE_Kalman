function return_val = Newton_method(start,df,d2f,tol,max_it)

pom = start;
for i = 1:max_it
    new_pom = pom - double(subs(df,pom))/double(subs(d2f,pom));
    if (abs(new_pom - pom) <= tol)
        return_val = new_pom;
    break;
    else
        pom = new_pom;
        return_val = new_pom;
    end
end

end