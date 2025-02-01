using System.Numerics;
using System.Diagnostics;

public static class Program
{
    private struct Ray(in Vector3 origin, in Vector3 direction)
    {
        public Vector3 Origin = origin;
        public Vector3 Direction = direction;
    }

    private enum Material
    {
        Diffuse = 0,
        Reflective = 1,
        Refractive = 2
    }

    private readonly struct Sphere(float radius, in Vector3 center, in Vector3 emission, in Vector3 color, Material material)
    {
        public Vector3 Center { get; } = center;
        public float Radius { get; } = radius;
        public Vector3 Emission { get; } = emission;
        public Material Material { get; } = material;
        public Vector3 Color { get; } = color;

        public float Intersect(in Ray ray)
        {
            const float eps = 1e-3f;
            var f = ray.Origin - Center;
            var a = ray.Direction.LengthSquared();
            var b = -Vector3.Dot(f, ray.Direction);
            var c = f.LengthSquared() - Radius * Radius;
            var delta = Radius * Radius - (f + b / a * ray.Direction).LengthSquared();
            if (delta < 0f) return 0f;
            var q = b + MathF.Sign(b) * MathF.Sqrt(a * delta);
            var t0 = c / q;
            if (t0 > eps) return t0;
            var t1 = q / a;
            return t1 > eps ? t1 : 0f;
        }
    }
    
    private static readonly Sphere[] Spheres =
    [
        new (1e4f, new Vector3(1e4f + 1f, 40.8f, 81.6f), Vector3.Zero, new Vector3(0.75f, 0.25f, 0.25f),
            Material.Diffuse),
        new (1e4f, new Vector3(-1e4f + 99f, 40.8f, 81.6f), Vector3.Zero, new Vector3(0.25f, 0.25f, 0.75f),
            Material.Diffuse),
        new (1e4f, new Vector3(50f, 40.8f, 1e4f), Vector3.Zero, new Vector3(0.75f, 0.75f, 0.75f),
            Material.Diffuse),
        new (1e4f, new Vector3(50f, 40.8f, -1e4f + 181f), Vector3.Zero, Vector3.Zero,
            Material.Diffuse),
        new (1e4f, new Vector3(50f, 1e4f, 81.6f), Vector3.Zero, new Vector3(0.75f, 0.75f, 0.75f),
            Material.Diffuse),
        new (1e4f, new Vector3(50f, -1e4f + 81.6f, 81.6f), Vector3.Zero, new Vector3(0.75f, 0.75f, 0.75f),
            Material.Diffuse),
        new (16.5f, new Vector3(27f, 16.5f, 47f), Vector3.Zero, new Vector3(0.999f, 0.999f, 0.999f),
            Material.Reflective),
        new (16.5f, new Vector3(73f, 16.5f, 78f), Vector3.Zero, new Vector3(0.999f, 0.999f, 0.999f),
            Material.Refractive),
        new (600f, new Vector3(50f, 681.6f - 0.27f, 81.6f), new Vector3(12f, 12f, 12f), Vector3.Zero,
            Material.Diffuse)
    ];

    private static bool Intersect(in Ray ray, ref float t, ref int id)
    {
        t = float.PositiveInfinity;
        for (var i = 0; i < Spheres.Length; ++i)
        {
            var d = Spheres[i].Intersect(ray);
            if (d == 0f || d >= t) continue;
            t = d;
            id = i;
        }

        return !float.IsPositiveInfinity(t);
    }

    private static Vector3 Radiance(Ray ray, Random rng)
    {
        const int maxDepth = 8;
        const float eps = 1e-2f;
        var t = 0f;
        var id = 0;
        var f = Vector3.One;
        var r = Vector3.Zero;
        for (var depth = 0; depth < maxDepth; ++depth)
        {
            if (!Intersect(ray, ref t, ref id)) break;
            var obj = Spheres[id];
            var x = ray.Origin + t * ray.Direction;
            var n = Vector3.Normalize(x - obj.Center);
            var nl = Vector3.Dot(n, ray.Direction) < 0f ? n : -n;
            r += f * obj.Emission;
            f *= obj.Color;
            var p = MathF.Max(f.X, MathF.Max(f.Y, f.Z));
            const int rrDepth = 5;
            if (depth > rrDepth)
            {
                if (rng.NextSingle() >= p)
                {
                    break;
                }

                f /= p;
            }

            switch (obj.Material)
            {
                case Material.Diffuse:
                {
                    var r1 = 2f * MathF.PI * rng.NextSingle();
                    var r2 = rng.NextSingle();
                    var r2Sqrt = MathF.Sqrt(r2);
                    var u = Vector3.Normalize(Vector3.Cross(MathF.Abs(nl.X) > 0.1f ? Vector3.UnitY : Vector3.UnitX,
                        nl));
                    var v = Vector3.Cross(nl, u);
                    ray.Origin = x + eps * nl;
                    ray.Direction = (u * MathF.Cos(r1) + v * MathF.Sin(r1)) * r2Sqrt + nl * MathF.Sqrt(1f - r2);
                    break;
                }
                case Material.Reflective:
                {
                    ray.Origin = x + eps * nl;
                    ray.Direction -= n * 2f * Vector3.Dot(n, ray.Direction);
                    break;
                }
                case Material.Refractive:
                default:
                {
                    var rDir = ray.Direction - n * 2f * Vector3.Dot(n, ray.Direction);
                    var into = Vector3.Dot(n, nl) > 0f;
                    const float nc = 1f;
                    const float nt = 1.5f;
                    var nnt = into ? nc / nt : nt / nc;
                    var ddn = Vector3.Dot(ray.Direction, nl);
                    var cos2T = 1f - nnt * nnt * (1f - ddn * ddn);
                    if (cos2T < 0f)
                    {
                        ray.Origin = x + eps * nl;
                        ray.Direction = rDir;
                        break;
                    }

                    var tDir = ray.Direction * nnt - n * ((into ? 1f : -1f) * (ddn * nnt + MathF.Sqrt(cos2T)));
                    const float a = nt - nc;
                    const float b = nt + nc;
                    const float r0 = a * a / (b * b);
                    var c = 1f - (into ? -ddn : Vector3.Dot(tDir, n));
                    var re = r0 + (1f - r0) * c * c * c * c * c;
                    var tr = 1f - re;
                    var rrThreshold = 0.25f + 0.5f * re;
                    if (rng.NextSingle() < rrThreshold)
                    {
                        ray.Origin = x + eps * nl;
                        ray.Direction = rDir;
                        f *= re / rrThreshold;
                    }
                    else
                    {
                        ray.Origin = x - eps * nl;
                        ray.Direction = tDir;
                        f *= tr / (1f - rrThreshold);
                    }

                    break;
                }
            }
        }

        return r;
    }

    private static int ToInt(float x) => Convert.ToInt32(MathF.Pow(x, 1f / 2.2f) * 255.0f + 0.5f);

    public static void Main(string[] args)
    {
        const int w = 1024, h = 768;
        var samples = args.Length > 0 ? int.Parse(args[0]) : 64;
        var c = new Vector3[w * h];
        Action<int, Random> renderRow = (y, rng) =>
        {
            var cam = new Ray(new Vector3(50f, 52f, 295.6f), Vector3.Normalize(new Vector3(0f, -0.042612f, -1f)));
            var cx = new Vector3(w * 0.5135f / h, 0f, 0f);
            var cy = Vector3.Normalize(Vector3.Cross(cx, cam.Direction)) * 0.5135f;
            for (var x = 0; x < w; ++x)
            {
                for (var sy = 0; sy < 2; ++sy)
                {
                    for (var sx = 0; sx < 2; ++sx)
                    {
                        var r = Vector3.Zero;
                        for (var s = 0; s < samples; ++s)
                        {
                            var r1 = 2f * rng.NextSingle();
                            var dx = r1 < 1f ? MathF.Sqrt(r1) - 1f : 1f - MathF.Sqrt(2f - r1);
                            var r2 = 2f * rng.NextSingle();
                            var dy = r2 < 1f ? MathF.Sqrt(r2) - 1f : 1f - MathF.Sqrt(2f - r2);
                            var d = cx * (((sx + 0.5f + dx) / 2f + x) / w - 0.5f) +
                                    cy * (((sy + 0.5f + dy) / 2f + y) / h - 0.5f) + cam.Direction;
                            r += Radiance(new Ray(cam.Origin + d * 140f, Vector3.Normalize(d)), rng) * (1f / samples);
                        }

                        c[x + (h - 1 - y) * w] += 0.25f * Vector3.Clamp(r, Vector3.Zero, Vector3.One);
                    }
                }
            }
        };
        var sw = Stopwatch.StartNew();
        Parallel.For(0, h, y =>
        {
            var rng = new ThreadLocal<Random>(() => new Random(y));
            Debug.Assert(rng.Value != null, "rng.Value != null");
            renderRow(y, rng.Value);
        });
        sw.Stop();
        Console.WriteLine($"Elapsed time: {sw.Elapsed}");
        using var f = File.CreateText("image.ppm");
        f.Write($"P3\n{w} {h}\n255\n");
        for (var i = 0; i < c.Length; ++i)
        {
            f.Write($"{ToInt(c[i].X)} {ToInt(c[i].Y)} {ToInt(c[i].Z)} ");
        }
    }
}