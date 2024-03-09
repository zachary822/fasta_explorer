/** @type {import('next').NextConfig} */
const nextConfig = {
  async rewrites() {
    return [
      {
        source: "/api/:path*",
        destination: `${process.env.API_URL ?? "http://api:8000"}/:path*`, // Proxy to Backend
      },
    ];
  },
};

export default nextConfig;
