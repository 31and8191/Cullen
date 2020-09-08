#include <pari/pari.h>

#define prec 100
#define BOUND 10
#define BOUND2 1000000
static GEN a;
static GEN t;
static GEN c;
static GEN z;
/*End of global vars*/

main()	  /* void */
{
  GEN p1 = gen_0;	  /* vec */
  pari_init(4000000000,100000);
  pari_sp ltop = avma;
  a = pol_x(fetch_user_var("a"));
  t = pol_x(fetch_user_var("t"));
  c = pol_x(fetch_user_var("c"));
  z = pol_x(fetch_user_var("z"));
  {
    long l2;
    p1 = cgetg(701, t_VEC);
    for (l2 = 1; l2 <= 700; ++l2)
      gel(p1, l2) = gen_0;
  }
  a = p1;
  {
    pari_sp btop = avma;
    long j;
    for (j = 2; j <= 700; ++j)
    {
      gel(a, j) = gdiv(gsubgs(gpowgs(t, j), 1), gsubgs(t, 1));
         if (gc_needed(btop, 1))
        a = gerepilecopy(btop, a);
     }
  }
  {
    pari_sp btop = avma;
    long s;
    for (s = 2; s <= BOUND2; ++s)
    {
      pari_printf("%ld\n", s);
      {
        pari_sp btop = avma;
        long n;
        for (n = 1; n <= BOUND; ++n)
        {
         {
          GEN p3 = gen_0;
          c = gaddgs(gmulsg(n, powis(stoi(s), n)), 1);
          p3 = gaddgs(gdiv(glog(c, prec), glog(gen_2, prec)), 1);
          {
            pari_sp btop = avma;
            GEN y = gen_0;
            for (y = stoi(3); gcmp(y, p3) <= 0; y = gaddgs(y, 1))
            {
              GEN p4;	  /* vec */
              int num;
              GEN f = gsub(gel(a, gtos(y)),c);
              GEN smallprime=stoi(17);
              while(gtos(smallprime)<100)
              {
                num = glength(polrootsmod(f,smallprime));
                smallprime=nextprime(gaddgs(smallprime,1));
                if (num==0)
                  smallprime=stoi(100);
              }
if(gtos(smallprime)!=100)
{
              p4 = cgetg(3, t_VEC);
              gel(p4, 1) = gen_0;
              gel(p4, 2) = gcopy(c);
              z = compo(realroots(gsub(gel(a, gtos(y)), c), p4, prec),1);
              if (!gequal0(gsub(gsubst(gel(a, gtos(y)), gvar(t), gfloor(z)), c)))
              {
                if (gequal0(gsub(gsubst(gel(a, gtos(y)), gvar(t), gceil(z)), c)))
                  pari_printf("s=%ld n=%ld y= %Ps\n", s, n, y);
              }
              else
                pari_printf("s=%ld n=%ld y= %Ps\n", s, n, y);
}
            }
                            if (gc_needed(btop, 1))
                  gerepileall(btop, 2, &y, &z);
          }
         }
                   if (gc_needed(btop, 1))
            gerepileall(btop, 2, &c, &z);
        }
      }
                        if (gc_needed(btop, 1))
            gerepileall(btop, 2, &c, &z);
    }
  }
    gerepileall(ltop, 3, &a, &c, &z);
    return 0;
}

