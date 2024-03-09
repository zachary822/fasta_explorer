import Upload from "@/components/Upload";

export default function Home() {
  return (
    <main className="bg-secondary">
      <div className="flex min-h-screen flex-col items-center md:container md:px-auto p-12">
        <Upload />
      </div>
    </main>
  );
}
