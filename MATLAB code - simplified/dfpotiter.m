function dfpotiter()
%c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%c this evaluates the two electron fock operator.
%c p(i,j) is the density matrix.
%c------------------------------------------------------------------
global nsym nbas nocorb ...
    nskipe nskipc ...
    df_two orbjc vec mbasis gl gs orbje rp

Jll(mbasis,mbasis)=0;
Jss(mbasis,mbasis)=0;
Kll(mbasis,mbasis)=0;
Kss(mbasis,mbasis)=0;
Kls(mbasis,mbasis)=0;
Ksl(mbasis,mbasis)=0;

for isym=1:nsym
    for ip=1:nbas(isym)

        jp=ip+nskipe(isym);
        kp=ip+2*nskipe(isym);
        indexp=ip+nbas(isym)+2*nskipe(isym);
        j1 = orbje(jp);

        for iq=1:nbas(isym)

            jq=iq+nskipe(isym);
            kq=iq+2*nskipe(isym);
            indexq=iq+nbas(isym)+2*nskipe(isym);
            j2 = orbje(jq);

            for jsym=1:nsym
                for ik=1:nocorb(jsym)

                    jk=ik+nskipc(jsym);
                    kk=ik+2*nskipe(jsym);

                    Nk = 2*orbjc(jk)+1;

                    for ir=1:nbas(jsym)

                        jr=ir+nskipe(jsym);
                        kr=ir+2*nskipe(jsym);
                        indexr=ir+nbas(jsym)+2*nskipe(jsym);

                        for is=1:nbas(jsym)

                            js=is+nskipe(jsym);
                            ks=is+2*nskipe(jsym);
                            indexs=is+nbas(jsym)+2*nskipe(jsym);

                            Dll = Nk*vec(kr,kk)*vec(ks,kk);
                            Dss = Nk*vec(indexr,kk)*vec(indexs,kk);
                            Dls = Nk*vec(kr,kk)*vec(indexs,kk);
                            Dsl = Nk*vec(indexr,kk)*vec(ks,kk);

                            Jllll = valint(gl(jp,:).*gl(jq,:).*valint2(0,gl(jr,:),gl(js,:)).*rp);
                            Jllss = valint(gl(jp,:).*gl(jq,:).*valint2(0,gs(jr,:),gs(js,:)).*rp);
                            Jssss = valint(gs(jp,:).*gs(jq,:).*valint2(0,gs(jr,:),gs(js,:)).*rp);
                            Jssll = valint(gs(jp,:).*gs(jq,:).*valint2(0,gl(jr,:),gl(js,:)).*rp);

                            Jll(kp,kq)=Jll(kp,kq)+Dll*Jllll+Dss*Jllss;
                            Jss(indexp,indexq)=Jss(indexp,indexq)+Dss*Jssss+Dll*Jssll;

                            for v=abs(j1-j2):(j1+j2)

                                factor = Wigner3j(j1,v,j2,0.5,0,-0.5)^2;
                                %Vladimir Sovkov (2022). Wigner 3j-6j-9j ...
                                %https://www.mathworks.com/matlabcentral/fileexchange/74069-wigner-3j-6j-9j

                                Kllll = valint(gl(jp,:).*gl(jr,:).*valint2(v,gl(jq,:),gl(js,:)).*rp);
                                Kssss = valint(gs(jp,:).*gs(jr,:).*valint2(v,gs(jq,:),gs(js,:)).*rp);
                                Klsls = valint(gl(jp,:).*gl(jr,:).*valint2(v,gs(jq,:),gs(js,:)).*rp);
                                Kslsl = valint(gs(jp,:).*gs(jr,:).*valint2(v,gl(jq,:),gl(js,:)).*rp);

                                Kll(kp,kq)=Kll(kp,kq)+factor*Dll*Kllll;
                                Kss(indexp,indexq)=Kss(indexp,indexq)+factor*Dss*Kssss;
                                Kls(kp,indexq)=Kls(kp,indexq)+factor*Dls*Klsls;
                                Ksl(indexp,kq)=Ksl(indexp,kq)+factor*Dsl*Kslsl;

                            end
                        end
                    end
                end
            end
        end
    end
end
J = Jll + Jss;
K = Kll + Kss + Ksl + Kls;
df_two = J-K;