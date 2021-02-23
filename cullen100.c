/*-*- compile-command: "/u01/grantham/gp2c/gp2c-0.0.11pl4/./config/missing "; -*-*/
#include <pari/pari.h>

#define prec 100
static GEN a;
static GEN t;
static GEN c;
static GEN z;
/*End of global vars*/

main(int argc, char *argv[])	  /* void */
{
  if(argc!=3)
  {
    printf("Need 2 arguments. Had %d.\n",argc-1);
    exit(99);
  }
  long n0,s;
  s=atoi(argv[1]);
  n0=atoi(argv[2]);
  GEN p1 = gen_0;	  /* vec */
  pari_init(400000000000,100000);
  pari_sp ltop = avma;
  a = pol_x(fetch_user_var("a"));
  t = pol_x(fetch_user_var("t"));
  c = pol_x(fetch_user_var("c"));
  z = pol_x(fetch_user_var("z"));
  {
    long l2;
    p1 = cgetg(10000, t_VEC);
    for (l2 = 1; l2 <= 10000; ++l2)
      gel(p1, l2) = gen_0;
  }
  a = p1;
  {
    long j;
    for (j = 2; j <= 10000; ++j)
    {
      gel(a, j) = gdiv(gsubgs(gpowgs(t, j), 1), gsubgs(t, 1));
     }
  }
  {
    {
      printf("s=%ld\n", s);
      {
        long n;
        for (n = 1; n <= n0; ++n)
        {
         printf("n=%d\n",n);
         fflush(stdout);
         {
          GEN p3 = gen_0;
          c = gaddgs(gmulsg(n, powis(stoi(s), n)), 1);
          p3 = gaddgs(gdiv(glog(c, prec), glog(gen_2, prec)), 1);
          //p3 = gmin(p3,stoi(s));
          {
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
                //if(!geq(gceil(z),stoi(1)))
                  pari_printf("s=%ld n=%ld q=%Ps b=%Ps\n", s, n, y,gceil(z));
              }
              else
                //if(!geq(gfloor(z),stoi(1)))
                pari_printf("s=%ld n=%ld q=%Ps b=%Ps\n", s, n, y, gfloor(z));
}
            }
          }
         }
        }
      }
    }
  }
    printf("Done.\n");
    return 0;
}

